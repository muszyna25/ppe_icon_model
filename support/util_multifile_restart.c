#include <assert.h>
#include <stdbool.h>
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

//A wrapper for opendir() which outputs an appropriate error message and returns NULL on failure.
DIR* openDirectory(const char* path) {
	errno = 0;
	DIR* result = opendir(path);
	if(result) return result;

	//the call failed, print the appropriate error message
	switch(errno) {
		case EACCES:
			fprintf(stderr, "error: could not read directory at \"%s\"\n", path);
			return NULL;

		case EMFILE:	//fallthrough
		case ENFILE:	//fallthrough
		case ENOMEM:
			fprintf(stderr, "error: insufficient resources to read directory: %s\n", strerror(errno));
			return NULL;

		case ENOENT:
			//If this condition is met, someone has actively deleted the directory after we created it,
			//and might be racing with us creating some symlink to let us overwrite some data of their choosing.
			//Stop this by failing noisily.
			fprintf(stderr, "fatal error: directory disappeared after creation. Is somebody trying to do something nasty with your process?\n");
			abort();

		case ENOTDIR:
			fprintf(stderr, "error: path at \"%s\" exists and is not a directory\n", path);
			return NULL;

		default:
			fprintf(stderr, "fatal error: unknown error reading directory at \"%s\": \"%s\"\n", path, strerror(errno));
			abort();
	}

	//unreachable code
	abort();
}

bool isPrefix(const char* prefix, const char* string) {
	for(; *prefix; prefix++, string++) {
		if(*prefix != *string) return false;
	}
	return true;
}

typedef enum {
	kStdFile,	//this is used for the . and .. entries in the directory
	kAttributesFile,
	kMetadataFile,
	kPayloadFile,
	kUnknownFile
} FileClass;

FileClass classifyFileName(const char* name, int* out_patchId, int* out_procId) {
	assert(out_patchId);
	assert(out_procId);

	*out_patchId = *out_procId = -1;
	if(!strcmp(".", name)) return kStdFile;
	if(!strcmp("..", name)) return kStdFile;
	if(!strcmp("attributes.nc", name)) return kAttributesFile;

	if(!isPrefix("patch", name)) return kUnknownFile;
	name += strlen("patch");	//skip the constant prefix

	//get the patch ID
	/*const*/ char* endPtr;	//the function signature of strtol() is broken, so we need to skip the const
	long patchId = strtol(name, &endPtr, 10);
	if(endPtr == name) return kUnknownFile;	//no patch ID given at the expected place
	if(patchId != (long)(int)patchId) return kUnknownFile;	//patch ID out of range of integer
	if(patchId < 1) return kUnknownFile;	//patch IDs start at 1
	*out_patchId = (int)patchId;
	name = endPtr;	//skip the patch ID

	if(!strcmp("_metadata", name)) return kMetadataFile;
	if(!isPrefix("_", name)) return kUnknownFile;
	name += strlen("_");	//skip the constant prefix

	//get process ID
	long procId = strtol(name, &endPtr, 10);
	if(endPtr == name) return kUnknownFile;	//no process ID given at the expected place
	if(procId != (long)(int)procId) return kUnknownFile;	//process ID out of range of integer
	if(procId < 0) return kUnknownFile;	//process IDs are positive
	*out_procId = (int)procId;
	name = endPtr;	//skip the process ID

	if(!strcmp(".nc", name)) return kPayloadFile;
	return kUnknownFile;
}

typedef struct FileParams {
	int patchId, procId;
} FileParams;

//allow sorting of FileParams arrays via qsort()
int compareFileParams(const void* aArg, const void* bArg) {
	const FileParams* a = aArg, *b = bArg;
	if(a->patchId != b->patchId) return a->patchId - b->patchId;
	if(a->procId != b->procId) return a->procId - b->procId;
	return 0;
}

//this may reorder the contents of the given arrays as a side effect
//
//if the expected domain or file count is unknown, just pass 0 in the respective argument
bool isFileListConsistent(size_t numof_metadataFiles, FileParams* metadataFiles, size_t numof_payloadFiles, FileParams* payloadFiles, int expectedDomainCount, int expectedFileCount) {
	//sort the two arrays, so we can easily compare them element by element
	qsort(metadataFiles, numof_metadataFiles, sizeof(*metadataFiles), compareFileParams);
	qsort(payloadFiles, numof_payloadFiles, sizeof(*payloadFiles), compareFileParams);

	//if we expect a certain domain count, check that the metadata files are consistent with it
	if(expectedDomainCount) {
		int expectedDomain = 1;
		for(int metadataIndex = 0; metadataIndex < numof_metadataFiles; metadataIndex++, expectedDomain++) {
			assert(metadataFiles[metadataIndex].patchId >= expectedDomain);	//we sorted the array, so this should be true
			if(metadataFiles[metadataIndex].patchId > expectedDomain) {
				fprintf(stderr, "warning: domain %d missing from restart multifile\n", expectedDomain);
				metadataIndex--;	//recheck this entry with the next expected domain
			}
		}
		int lastDomain = expectedDomain - 1;
		if(lastDomain < expectedDomainCount) fprintf(stderr, "warning: domains %d through %d missing from restart multifile\n", lastDomain + 1, expectedDomainCount);
		if(lastDomain > expectedDomainCount) {
			fprintf(stderr, "error: restart multifile contains unexpected domain (found metadata file for domain %d, but expected only %d domains)\n", lastDomain, expectedDomainCount);
			return false;
		}
	}

	//if we expect a certain file count, check that the total count of payload files is as expected
	if(expectedFileCount) {
		if(numof_metadataFiles*expectedFileCount != numof_payloadFiles) {
			fprintf(stderr, "error: unexpected number of payload files in restart multifile (expected %d files, %d files per patch, but found %d payload files)\n", numof_metadataFiles*expectedFileCount, expectedFileCount, numof_payloadFiles);
			return false;
		}
	}

	//walk through the payload files and check that they have a corresponding metadata file
	for(int metadataIndex = 0, payloadIndex = 0; payloadIndex < numof_payloadFiles; payloadIndex++) {
		//advance the metadataIndex so that it points to the respective metadata file
		while(metadataIndex < numof_metadataFiles && metadataFiles[metadataIndex].patchId < payloadFiles[payloadIndex].patchId) metadataIndex++;
		if(metadataIndex == numof_metadataFiles) return false;	//the list of metadata files was exhausted before the patch ID of the payload file was reached

		if(metadataFiles[metadataIndex].patchId != payloadFiles[payloadIndex].patchId) return false;	//there is no matching metadata file for the payload file

		//check the expected range of the procId
		if(payloadFiles[payloadIndex].procId < 0) return false;	//we never expect files with a negative procId
		if(expectedFileCount && payloadFiles[payloadIndex].procId >= expectedFileCount) return false;	//we can only check the upper bound if we are expecting a certain file count
	}
	return true;
}

//if the expected domain or file count is unknown, just pass 0 in the respective argument
int checkDirContents(DIR* directory, const char* dirPath, int expectedDomainCount, int expectedFileCount) {
	bool haveAttributesFile = false;
	size_t numof_metadataFiles = 0, sizeof_metadataFiles = 8;
	FileParams* metadataFiles = malloc(sizeof_metadataFiles*sizeof(*metadataFiles));
	size_t numof_payloadFiles = 0, sizeof_payloadFiles = 8;
	FileParams* payloadFiles = malloc(sizeof_payloadFiles*sizeof(*payloadFiles));

	if(!metadataFiles || !payloadFiles) fprintf(stderr, "memory allocation failure\n"), abort();

	//get the underlying file descriptor of the directory, we need this so we can use fstatat() rather than lstat(), which should provide a significant performance benefit
	int dirFileDescriptor = dirfd(directory);
	if(dirFileDescriptor < 0) fprintf(stderr, "fatal error: dirfd() returned an error: \"%s\"", strerror(errno)), abort();

	struct dirent* curEntry;
	while((curEntry = readdir(directory))) {
		int patchId, procId;
		FileClass type = classifyFileName(curEntry->d_name, &patchId, &procId);
		switch(type) {
			case kStdFile: break;	//ignore the . and .. entries

			case kAttributesFile:
				assert(!haveAttributesFile);
				haveAttributesFile = true;
				break;

			case kMetadataFile:
				if(numof_metadataFiles == sizeof_metadataFiles) {
					metadataFiles = realloc(metadataFiles, (sizeof_metadataFiles *= 2)*sizeof(*metadataFiles));
					if(!metadataFiles) fprintf(stderr, "memory allocation failure\n"), abort();
				}
				metadataFiles[numof_metadataFiles++] = (FileParams){
					.patchId = patchId,
					.procId = -1
				};
				break;

			case kPayloadFile:
				if(numof_payloadFiles == sizeof_payloadFiles) {
					payloadFiles = realloc(payloadFiles, (sizeof_payloadFiles *= 2)*sizeof(*payloadFiles));
					if(!payloadFiles) fprintf(stderr, "memory allocation failure\n"), abort();
				}
				payloadFiles[numof_payloadFiles++] = (FileParams){
					.patchId = patchId,
					.procId = procId
				};
				break;

			case kUnknownFile:
				fprintf(stderr, "error: directory at \"%s\" is not a valid multifile: unknown file present \"%s\"\n", dirPath, curEntry->d_name);
				return -1;
		}
		if(type != kStdFile) {
			//check that the file is a regular file
			struct stat fileInfo;
			if(fstatat(dirFileDescriptor, curEntry->d_name, &fileInfo, AT_SYMLINK_NOFOLLOW)) {
				fprintf(stderr, "error stating file \"%s\" within directory \"%s\"\n", curEntry->d_name, dirPath);
				return -1;
			}
			if((fileInfo.st_mode & S_IFMT) != S_IFREG) {
				fprintf(stderr, "error: directory at \"%s\" is not a valid multifile: \"%s\" is not a regular file\n", dirPath, curEntry->d_name);
				return -1;
			}
		}
	}
	if(!haveAttributesFile) {
		fprintf(stderr, "error: directory at \"%s\" is not a valid multifile: no \"attributes.nc\" found\n", dirPath);
		return -1;
	}
	if(!isFileListConsistent(numof_metadataFiles, metadataFiles, numof_payloadFiles, payloadFiles, expectedDomainCount, expectedFileCount)) {
		//XXX: I'm not sure what to do in this case: There are no unknown files present, but the directory is also not a valid multifile.
		//     It could be a corrupted one (likely), or it could be a directory that just happens to contain a file `attributes.nc` and one or more `patchN_M.nc` files (unlikely).
		//     In the former case, the right course of action would be to throw a warning that we are overwriting a corrupted multifile, and proceed normally,
		//     in the later case, the right course of action is to abort the operation to avoid overwriting data that we don't know about.
		//     Either choice is wrong as it can lead to data loss (either by overwriting data we didn't produce, or by aborting a run we needn't abort).
		//
		//     Since the corrupted multifile explanation seems much more likely to me, I opt for the warn-and-overwrite approach here.
		fprintf(stderr, "warning: inconsistent multifile detected at \"%s\"\n", dirPath);
	}

	//cleanup
	free(metadataFiles);
	free(payloadFiles);
	return 0;
}

int checkMultifileDir(const char* path, int expectedDomainCount, int expectedFileCount) {
	DIR* directory = openDirectory(path);
	if(!directory) return -1;
	int result = checkDirContents(directory, path, expectedDomainCount, expectedFileCount);
	if(closedir(directory)) {
		fprintf(stderr, "fatal error: can't close directory, so something must be very wrong\n");
		abort();
	}
	return result;
}

int cleanDirContents(DIR* directory, const char* dirPath) {
	//get the underlying file descriptor of the directory, we need this so we can use unlinkat() rather than unlink(), which should provide a significant performance benefit
	int dirFileDescriptor = dirfd(directory);
	if(dirFileDescriptor < 0) fprintf(stderr, "fatal error: dirfd() returned an error: \"%s\"", strerror(errno)), abort();

	struct dirent* curEntry;
	while((curEntry = readdir(directory))) {
		int patchId, procId;
		FileClass type = classifyFileName(curEntry->d_name, &patchId, &procId);
		switch(type) {
			case kStdFile: break;	//ignore the . and .. entries

			case kAttributesFile:	//fallthrough
			case kMetadataFile:	//fallthrough
			case kPayloadFile:
				if(unlinkat(dirFileDescriptor, curEntry->d_name, 0)) {
					fprintf(stderr, "fatal error: could not unlink file \"%s\" from directory \"%s\" (%s)\n", curEntry->d_name, dirPath, strerror(errno));
					return -1;
				}
				break;

			case kUnknownFile:
				fprintf(stderr, "fatal error: unknown file present (file \"%s\" in directory \"%s\") that was not detected in the previous check. This means we were either racing with some other process creating the file, or there is a bug in the sanity checking code.\n", curEntry->d_name, dirPath);
				return -1;

			default:
				fprintf(stderr, "assertion failed: unreachable code reached\n");
				abort();
		}
	}
	return 0;
}

int createEmptyMultifileDir(const char* path) {
	errno = 0;
	if(mkdir(path, 0777)) {
		//ok, so we couldn't create the directory, but for which reason?
		switch(errno) {
			case EEXIST: {
				//This is not a fatal case: it may be the case that we are overwriting an existing restart multifile.
				//So we just check that we really have a directory, and that it contains only files that would belong to a restart multifile.
				//If that turns out to be the case, we clear the contents of the directory before continuing.
				DIR* directory = openDirectory(path);
				if(!directory) return -1;
				if(checkDirContents(directory, path, 0, 0)) return -1;
				rewinddir(directory);
				if(cleanDirContents(directory, path)) return -1;
				if(closedir(directory)) {
					fprintf(stderr, "fatal error: can't close directory, so something must be very wrong\n");
					abort();
				}
				return 0;
			}

			case EACCES:	//fallthrough
			case EPERM:	//fallthrough
			case EROFS:
				fprintf(stderr, "error: insufficient permissions for creating directory at \"%s\", check that the directory is writeable\n", path);
				return -1;

			case EFAULT:
				fprintf(stderr, "assertion failed: invalid pointer passed to createEmptyMultifile()\n");
				abort();

			case ELOOP:
				fprintf(stderr, "error: too many symlinks encountered while trying to create directory at \"%s\"\n", path);
				return -1;

			case ENAMETOOLONG:
				fprintf(stderr, "error: path for restart multifile is too long (%d characters: \"%s\")\n", strlen(path), path);
				return -1;

			case ENOENT:
				fprintf(stderr, "error: a path component in \"%s\" does not exist or is a dangling symlink\n", path);
				return -1;

			case ENOMEM:
				fprintf(stderr, "fatal error: kernel out of memory\n");
				abort();

			case ENOSPC:
				fprintf(stderr, "error: could not create directory for restart multifile: out of storage space\n");
				return -1;

			case ENOTDIR:
				fprintf(stderr, "error: a path component in \"%s\" is not a directory\n", path);
				return -1;

			default:
				fprintf(stderr, "fatal error: unknown error while creating directory for restart multifile: \"%s\"\n", strerror(errno));
				abort();
		}
	}
	return 0;
}

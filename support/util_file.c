#define XOPEN_SOURCE 500

#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#define PID_LEN 8

int util_islink(char *path)
{
  char *buf;
  ssize_t len;
  long max_buf_len;
  int iret;

  iret = -1;

  max_buf_len = pathconf("/", _PC_NAME_MAX);

  buf = (char *) malloc(max_buf_len*sizeof(char));

  if ((len = readlink(path, buf, max_buf_len-1)) != -1)
    {
      buf[len] = '\0';
      /* file is a link ..., return C true */
      iret = 1;
    }
  else
    {
      if (errno == EINVAL)
        {
          /* file is not a link ... */
          iret = 0;
        }
    }

  free(buf);

  /* something else went wrong ... */
  return iret;
}

int util_tmpnam_len(void)
{
  return (L_tmpnam + 1 + PID_LEN);
}

int util_tmpnam(char *filename)
{
  char *ptr;
  char pid_string[PID_LEN + 1];
  pid_t pid;

  pid = getpid();

  sprintf(pid_string, "%8.8ld", (long) pid);

  ptr = tmpnam(filename);

  strcat(filename, pid_string);

  return strlen(filename);
}

long int util_filesize(char *filename)
{
  struct stat statbuf;

  if (stat(filename, &statbuf) == -1)
    {
      return 0;
    }

  return ((long int) statbuf.st_size);
}

int util_file_is_writable(char *filename)
{
  int result = 0;
  int rval = access (filename, W_OK);
  if (rval == 0)
    return 1;
  else if ((errno == EACCES) ||  (errno == EROFS))
    return 0;
  return 0;
}

int putFile(const char* path, size_t dataSize, const char* data, int fileMode)
{
  errno = 0;
  int file = open(path, O_CREAT | O_TRUNC | O_WRONLY, fileMode);
  if(file < 0)
    {
      fprintf(stderr, "error opening file at \"%s\" for writing: \"%s\"\n", path, strerror(errno));
      return -1;
    }

  while(dataSize)
    {
      errno = 0;
      ssize_t writtenBytes = write(file, data, dataSize);
      if(writtenBytes > 0)
        {
          dataSize -= writtenBytes;
          data += writtenBytes;
        }
      else
        {
          switch(errno)
            {
              case 0:	//fallthrough
              case EAGAIN:	//fallthrough
              case EINTR:
                  break;	//not fatal, just incomplete. Retry.

              default:
                  fprintf(stderr, "error while writing file at \"%s\": error code %d \"%s\"\n", path, errno, strerror(errno));
                  return -1;
            }
        }
    }

  return close(file);
}

//This wrapper to symlink() does not fail if there is already a symlink at linkName; it will simply overwrite the existing symlink.
//However, if something exists at linkName which is not a symlink, the error EEXIST is returned.
//
//May return EEXIST or any error returned by lstat(), unlink(), or symlink(), returns zero on success.
int createSymlink(const char* targetPath, const char* linkName)
{
  errno = 0;
  struct stat fileInfo;
  if(!lstat(linkName, &fileInfo))
    {
      if((fileInfo.st_mode & S_IFMT) == S_IFLNK)
        {
          if(unlink(linkName)) return errno;
        }
      else
        {
          return EEXIST;	//Something exists at linkName, and it's not a symlink. Bail out.
        }
    }
  else if(errno != ENOENT)	//ENOENT is the only non-fatal error, it just means that there is no such link there yet
    {
      return errno;
    }

  //At this point we know that there is no file at linkName.
  errno = 0;
  symlink(targetPath, linkName);
  return errno;
}


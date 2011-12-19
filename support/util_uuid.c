/*
 * This product is derived out of software developed by 
 *
 * The Apache Software Foundation (http://www.apache.org/).
 * 
 * Copyright 2011 Luis Kornblueh
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#define _POSIX_C_SOURCE 200112L

#ifdef __APPLE__
/* allow accessing syscall() */
#define _DARWIN_C_SOURCE
/* this prevents using the Darwin provided version of uuid */
#define _UUID_T 
#endif

#ifdef __linux__
#define _SVID_SOURCE	
#endif

/* required here for clean load of Darwin header files */
typedef struct {
  unsigned char data[16]; 
} uuid_t;

/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/time.h>

#ifndef _SX
#ifdef _AIX
#include <sys/thread.h>
#else
#include <sys/syscall.h>
#include <sys/types.h>
#endif
#endif

/****************************************************************************/

static unsigned short jrand_seed[3];

int get_random_fd(void)
{
  struct timeval  tv;
  static int      fd = -2;
  int             i;

  if (fd == -2) 
    {
      gettimeofday(&tv, 0);
      fd = open("/dev/urandom", O_RDONLY);
      if (fd == -1)
	fd = open("/dev/random", O_RDONLY | O_NONBLOCK);
      if (fd >= 0) {
	i = fcntl(fd, F_GETFD);
	if (i >= 0)
	  fcntl(fd, F_SETFD, i | FD_CLOEXEC);
      }
      srand((getpid() << 16) ^ getuid() ^ tv.tv_sec ^ tv.tv_usec);
      jrand_seed[0] = getpid() ^ (tv.tv_sec & 0xFFFF);
      jrand_seed[1] = getppid() ^ (tv.tv_usec & 0xFFFF);
      jrand_seed[2] = (tv.tv_sec ^ tv.tv_usec) >> 16;
    }

  gettimeofday(&tv, 0);
  for (i = (tv.tv_sec ^ tv.tv_usec) & 0x1F; i > 0; i--)
    rand();
  
  return fd;
}

void generate_random_bytes(void *buf, int nbytes)
{
  int i, n = nbytes, fd = get_random_fd();
  int lose_counter = 0;
  unsigned char *cp = (unsigned char *) buf;
  unsigned short tmp_seed[3];
  
  if (fd >= 0) 
    {
      while (n > 0) 
	{
	  i = read(fd, cp, n);
	  if (i <= 0) 
	    {
	      if (lose_counter++ > 16)
		break;
	      continue;
	    }
	  n -= i;
	  cp += i;
	  lose_counter = 0;
	}
    }

  for (cp = buf, i = 0; i < nbytes; i++)
    *cp++ ^= (rand() >> 7) & 0xFF;

  memcpy(tmp_seed, jrand_seed, sizeof(tmp_seed));

#if defined _AIX
  jrand_seed[2] = jrand_seed[2] ^ thread_self();
#elif defined _SX
  jrand_seed[2] = jrand_seed[2] ^ getpid();
#else
  jrand_seed[2] = jrand_seed[2] ^ syscall(SYS_gettid);
#endif
  
  for (cp = buf, i = 0; i < nbytes; i++)
    *cp++ ^= (jrand48(tmp_seed) >> 7) & 0xFF;

  memcpy(jrand_seed, tmp_seed, sizeof(jrand_seed)-sizeof(unsigned short));
}

#ifdef DEBUG

int main(int argc, char *argv[])
{
  unsigned char buffer[36];
  int i;

  generate_random_bytes(buffer, 36);

  for (i = 0; i < 36; i++)
    printf("%02x ", buffer[i]);
  printf("\n");

  return 0;
}

#endif

/****************************************************************************/

#define APR_USEC_PER_SEC APR_TIME_C(1000000)

uint64_t time_now(void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec * 1000000 + tv.tv_usec;
}

void get_system_time(uint64_t *uuid_time)
{
  *uuid_time = time_now();
  *uuid_time = (*uuid_time * 10) + (uint64_t) 0x01B21DD213814000;
}

void get_current_time(uint64_t *timestamp)
{
  uint64_t time_now;
  static uint64_t time_last = 0;
  static uint64_t fudge = 0;

  get_system_time(&time_now);
        
  if (time_last != time_now) 
    {
      if (time_last + fudge > time_now)
	fudge = time_last + fudge - time_now + 1;
      else
	fudge = 0;
      time_last = time_now;
    }
  else 
    {
      ++fudge;
    }
  
  *timestamp = time_now + fudge;
}

/****************************************************************************/

#define NODE_LENGTH 6

static int uuid_state_seqnum;
static unsigned char uuid_state_node[NODE_LENGTH] = { 0 };

int true_random(void)
{
  unsigned char buf[2];

  generate_random_bytes(buf, 2);
  return (buf[0] << 8) | buf[1];
}

void get_random_info(unsigned char node[NODE_LENGTH])
{
  generate_random_bytes(node, NODE_LENGTH);
}

void get_pseudo_node_identifier(unsigned char *node)
{
  get_random_info(node);
  node[0] |= 0x01;
}

void init_state(void)
{
    uuid_state_seqnum = true_random();
    get_pseudo_node_identifier(uuid_state_node);
}

void uuid_get(uuid_t *uuid)
{
  uint64_t timestamp;
  unsigned char *d = uuid->data;
  
  if (!uuid_state_node[0])
    init_state();
  
  get_current_time(&timestamp);
  
  /* time_low, uint32 */
  d[3] = (unsigned char)timestamp;
  d[2] = (unsigned char)(timestamp >> 8);
  d[1] = (unsigned char)(timestamp >> 16);
  d[0] = (unsigned char)(timestamp >> 24);
  /* time_mid, uint16 */
  d[5] = (unsigned char)(timestamp >> 32);
  d[4] = (unsigned char)(timestamp >> 40);
  /* time_hi_and_version, uint16 */
  d[7] = (unsigned char)(timestamp >> 48);
  d[6] = (unsigned char)(((timestamp >> 56) & 0x0F) | 0x10);
  /* clock_seq_hi_and_reserved, uint8 */
  d[8] = (unsigned char)(((uuid_state_seqnum >> 8) & 0x3F) | 0x80);
  /* clock_seq_low, uint8 */
  d[9] = (unsigned char)uuid_state_seqnum;
  /* node, byte[6] */
  memcpy(&d[10], uuid_state_node, NODE_LENGTH);
}

#define UUID_FORMATTED_LENGTH 36

void uuid_format(char *buffer, const uuid_t *uuid)
{
  const unsigned char *d = uuid->data;
  
  sprintf(buffer,
	  "%02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x",
	  d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7],
	  d[8], d[9], d[10], d[11], d[12], d[13], d[14], d[15]);
}

unsigned char parse_hexpair(const char *s)
{
  int result;
  
  if (isdigit(*s))
    {
      result = (*s - '0') << 4;
    }
  else
    {
      if (isupper(*s))
	{
	  result = (*s - 'A' + 10) << 4;
	}
      else
	{
	  result = (*s - 'a' + 10) << 4;
	}
    }
  
  ++s;
  if (isdigit(*s))
    {
      result |= (*s - '0');
    }
  else
    {
      if (isupper(*s))
	{
	  result |= (*s - 'A' + 10);
        }
      else
	{
	  result |= (*s - 'a' + 10);
        }
    }
  
  return (unsigned char) result;
}

int uuid_parse(uuid_t *uuid, const char *uuid_str)
{
  int i;
  unsigned char *d = uuid->data;

  for (i = 0; i < 36; ++i) 
    {
      char c = uuid_str[i];
      if (!isxdigit(c) && !(c == '-' && (i == 8 || i == 13 || i == 18 || i == 23)))
	return -1;
    }
  if (uuid_str[36] != '\0') 
    {
      return -1;
    }
  
  d[0] = parse_hexpair(&uuid_str[0]);
  d[1] = parse_hexpair(&uuid_str[2]);
  d[2] = parse_hexpair(&uuid_str[4]);
  d[3] = parse_hexpair(&uuid_str[6]);
  
  d[4] = parse_hexpair(&uuid_str[9]);
  d[5] = parse_hexpair(&uuid_str[11]);
  
  d[6] = parse_hexpair(&uuid_str[14]);
  d[7] = parse_hexpair(&uuid_str[16]);
  
  d[8] = parse_hexpair(&uuid_str[19]);
  d[9] = parse_hexpair(&uuid_str[21]);
  
  for (i = 6; i--;)
    d[10 + i] = parse_hexpair(&uuid_str[i*2+24]);

  return 0;
}

#ifdef DEBUG

int main(int argc, char *argv[])
{
  uuid_t uuid;
  uuid_t uuid2;
  char buf[UUID_FORMATTED_LENGTH + 1];

  uuid_get(&uuid);
  uuid_format(buf, &uuid);

  printf("%s \n", buf);

  uuid_parse(&uuid2, buf);
  if (memcmp(&uuid, &uuid2, sizeof(uuid)) != 0)
    {
      fprintf(stderr,"uuid_parse produced a different UUID ...\n");
    }

  uuid_get(&uuid);
  uuid_get(&uuid2);
  
  if (memcmp(&uuid, &uuid2, sizeof(uuid)) == 0)
    {
      fprintf(stderr,"uuid_get generated the same UUID twice ...\n");
    }

  return 0;
}

#endif

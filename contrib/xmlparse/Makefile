all:
	~/software/ragel-6.8/ragel/ragel xml_scan.rl
	gcc xml_scan.c -o xml_scan -std=gnu99 -march=native -O0 -g -Wall

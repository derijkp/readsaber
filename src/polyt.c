/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 */

#define _FILE_OFFSET_BITS 64
#define MAX_READNAME_SIZE 16384

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

int main(int argc, char *argv[]) {
	char readname[MAX_READNAME_SIZE];
	register int polyTstart,polyTend,polyTnum,polyTgaps,polyTpgaps,polyTpnum;
	register int polyAstart,polyAend,polyAnum,polyAgaps,polyApgaps,polyApnum;
	register int curpos;
	register char c;
	if ((argc != 1)) {
		fprintf(stderr,"Format is: polyt\n");
		exit(EXIT_FAILURE);
	}
	fprintf(stdout,"rname\trstart\trend\tsize\tstrand\tchromosome\n");
	while (fgets_unlocked(readname, MAX_READNAME_SIZE, stdin)) {
		/* remove newline */
		register char *crn=readname;
		while (*crn != '\0') {
			if (*crn == '\n' || *crn == '\t' || *crn == ' ') break;
			crn++;
		}
		*crn = '\0';
/* fprintf(stderr,"%s\n",readname); */
		curpos = 0;
		polyTstart = -1 ; polyTend = -1 ; polyTnum = 0 ; polyTgaps = 0 ; polyTpgaps = 0 ; polyTpnum = 0; 
		polyAstart = -1 ; polyAend = -1 ; polyAnum = 0 ; polyAgaps = 0 ; polyApgaps = 0 ; polyApnum = 0; 
		/* greedy approach, just streaming the sequence; */
		/* start at any T ; break of after > 3 diff chars, or 3 successive diff ; only output if >= 8 Ts detected */
		/* this will not allways get optimal polyA/T because it simply restarts from where broken of */
		while (1) {
			c=getc_unlocked(stdin);
/* fprintf(stderr,"%c\tcurpos=%d\tstart=%d\tnum=%d\tgaps=%d\tpnum=%d\tpgaps=%d\n",c,curpos,polyTstart,polyTnum,polyTgaps,polyTpnum,polyTpgaps); */
			if (c == 'T') {
				if (polyTstart != -1) {
					polyTnum++;
					/* reset polyTgaps to 1 if more than 5 successive Ts */
					if (polyTpnum >= 5 && polyTpgaps > 1) {polyTgaps = 1;}
				} else {
					polyTstart = curpos;
					polyTnum = 1;
					polyTgaps = 0;
					polyTpnum = 1;
				}
				polyTend = curpos+1;
				polyTpnum++;
				polyTpgaps = 0;
			} else if (polyTstart != -1) {
				polyTgaps++;
				if (polyTgaps > 3 || polyTpgaps >= 3 || c == '\n' || c == EOF) {
					if (polyTnum >= 8) {
						fprintf(stdout,"%s\t%d\t%d\t%d\t+\tpolyT\n",readname,polyTstart,polyTend,polyTend-polyTstart);
					}
					polyTstart = -1;
				}
				polyTpnum = 0;
				polyTpgaps++;
			}
			if (c == 'A') {
				if (polyAstart != -1) {
					polyAnum++;
					/* reset polyAgaps to 1 if more than 5 successive As */
					if (polyApnum >= 5 && polyApgaps > 1) {polyAgaps = 1;}
				} else {
					polyAstart = curpos;
					polyAnum = 1;
					polyAgaps = 0;
					polyApnum = 1;
				}
				polyAend = curpos+1;
				polyApnum++;
				polyApgaps = 0;
			} else if (polyAstart != -1) {
				polyAgaps++;
				if (polyAgaps > 3 || polyApgaps >= 3 || c == '\n' || c == EOF) {
					if (polyAnum >= 8) {
						fprintf(stdout,"%s\t%d\t%d\t%d\t-\tpolyT\n",readname,polyAstart,polyAend,polyAend-polyAstart);
					}
					polyAstart = -1;
				}
				polyApnum = 0;
				polyApgaps++;
			}
			if (c == '\n' || c == EOF) break;
			curpos++;
		}
		/* read (And ignore) the + and quality line */
		while (1) {
			c=getc_unlocked(stdin);
			if (c == '\n' || c == EOF) break;
		}
		while (1) {
			c=getc_unlocked(stdin);
			if (c == '\n' || c == EOF) break;
		}
	}
	exit(EXIT_SUCCESS);
}

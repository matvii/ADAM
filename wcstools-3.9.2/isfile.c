/* File isfile.c
 * July 9, 2014
 * By Jessica Mink, Harvard-Smithsonian Center for Astrophysics
 * Send bug reports to jmink@cfa.harvard.edu

   Copyright (C) 2014
   Smithsonian Astrophysical Observatory, Cambridge, MA USA

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 *
 * Return 1 if argument is a file, 2 if it is a directory, else 0
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libwcs/fitsfile.h"

static char *RevMsg = "ISFILE WCSTools 3.9.2, 15 May 2015, Jessica Mink (jmink@cfa.harvard.edu)";

int
main (ac, av)
int ac;
char **av;
{
    char *str;

    /* Check for version or help command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help")) {
	fprintf (stderr,"%s\n",RevMsg);
	fprintf (stderr,"Usage: Return 1 if argument is a file, 2 if a directory, else 0\n");
	exit (1);
	}
    else if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	fprintf (stderr,"%s\n",RevMsg);
	exit (1);
	}

    /* check to see if this is a file */
    else
	printf ("%d\n", isfile (str));

    exit (0);
}
/* Jul  9 2014	New program
 */

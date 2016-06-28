/*** simpos.c - search object by its name from command line arguments
 *** May 15, 2015
 *** By Jessica Mink, sort of after IPAC byname.c for searching NED
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "libwcs/wcs.h"
#include "libwcs/fitsfile.h"
#include "libwcs/wcscat.h"

extern int   ned_errno;
static char searchorder[4];
static int printall = 0;
static void PrintUsage();
static char *RevMsg = "SIMPOS 3.9.2, 15 May 2015, Jessica Mink SAO";

int
main (ac, av)
int  ac;
char *av[];
{
   
    int lobj;
    int i;
    int verbose = 0;
    int printid = 0;
    int printname = 0;
    int printdeg = 0;
    int nid = 0;
    int nobj = 0;
    int nfobj = 0;
    int cds = 0;
    int ned = 0;
    int viz = 0;
    int simbad = 0;
    int outsys = WCS_J2000;
    int sysj = WCS_J2000;
    int tabout = 0;
    double ra, dec;
    char *str, *objname, *posdec, *posra, *iend, *ieq;
    char *listfile;
    char rastr[32], decstr[32];
    char newobj[32];
    char *buff, *buffid, *idline, *posline, *errline, *id, *errend;
    char url[256];
    char *xbuff;
    int lbuff;
    FILE *flist;
    char cr, lf, eq;

    listfile = NULL;
    cr = (char) 13;
    lf = (char) 10;
    eq = '=';
    searchorder[0] = (char) 0;
    searchorder[1] = (char) 0;
    searchorder[2] = (char) 0;
    searchorder[3] = (char) 0;

    if ( !strcmp (*av, "cdspos") )
	cds = 1;
    if ( !strcmp (*av, "nedpos") ) {
	searchorder[0]='N';
	ned = 1;
	}
    if ( !strcmp (*av, "vizpos") ) {
	searchorder[0]='V';
	viz = 1;
	}

    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	PrintUsage (NULL);
    if (!strcmp (str, "version") || !strcmp (str, "-version"))
	PrintUsage ("version");


    /* crack arguments */
    for (av++; --ac > 0 && (**av == '-' || **av == '@'); av++) {

	if (**av == '-') {
            char c;
	    str = *av;
            while (c = *++str) {
            switch (c) {

        	case 'a':       /* Print all ID information */
	            printall = 1;
		    break;

        	case 'b':       /* Print coordinates in B1950 */
	            outsys = WCS_B1950;
		    break;

        	case 'c':       /* Use CDS server in France */
	            cds = 1;
		    break;

		case 'd':       /* Print coordinates in degrees */
		    printdeg++;
		    break;

		case 'e':       /* Print coordinates in ecliptic coordinates */
		    outsys = WCS_ECLIPTIC;
		    break;

		case 'g':       /* Print coordinates in galactic coordinates */
		    outsys = WCS_GALACTIC;
		    break;

		case 'i':       /* Print all IDs found by SIMBAD */
		    printid++;
		    break;

		case 'n':       /* Search NED database */
		    ned++;
		    strcat (searchorder, "N");
		    break;

		case 'o':       /* Print object name */
		    printname++;
		    break;

		case 's':       /* Database search order */
		    simbad++;
		    strcat (searchorder, "S");
		    break;

		case 't':       /* Print output in tab-separated table */
		    tabout++;
		    break;

		case 'v':       /* More verbosity, including first ID */
		    verbose++;
		    break;

		case 'z':       /* Search Vizier catalogs */
		    viz++;
		    strcat (searchorder, "V");
		    break;

		default:
		    PrintUsage(NULL);
		    break;
		}
		}
	    }

	/* File containing a list of object names */
	else if (**av == '@') {
	    listfile = *av + 1;
	    if (isfile (listfile)) {
		nfobj = getfilelines (listfile);
		if (verbose)
		    fprintf (stderr,"SIMPOS: %d objects from file %s\n",
			     nfobj, listfile);
		if ((flist = fopen (listfile, "r")) == NULL) {
		    fprintf (stderr,"SIMPOS: List file %s cannot be read\n",
			     listfile);
		    }
		}
	    else {
		printf ("SIMPOS: List file %s does not exist\n", listfile);
		listfile = NULL;
		}
	    
	    }
	}

    /* There are ac remaining file names starting at av[0] */
    if (ac == 0 && !listfile)
        PrintUsage (NULL);

    /* Set searches to SNV if not otherwise set */
    if ( searchorder[0] == (char) 0 )
	strcpy (searchorder, "SNV");
    if (listfile) {
	ac = nfobj;
	objname = newobj;
	}
    while (ac > 0) {

	if (listfile)
	    fgets (newobj, 32, flist);
	else
	    objname = *av++;
	ac--;

	/* Replace underscores and spaces with plusses */
	lobj = strlen (objname);
	for (i = 0; i < lobj; i++) {
	    if (objname[i] == '_')
		objname[i] = '+';
	    if (objname[i] == ' ') {
		if (i == lobj-1)
		    objname = (char) 0;
		else
		    objname[i] = '+';
		}
	    }
	if (verbose)
	    printf ("%s -> ", objname);

	for (i = 0; i < 256; i++) {
	    url[i] = (char) 0;
	    }
	if (cds)
	    strcpy (url, "http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame/-oI/");
	else
	    strcpy (url, "http://vizier.cfa.harvard.edu/viz-bin/nph-sesame/-oI/");
	strcat (url, searchorder);
	if ((printid || printall) && !strcmp (searchorder, "V"))
	    strcat (url, "A");
	strcat (url, "?");
	strcat (url, objname);
	if (verbose)
	    printf ("%s\n", url);
	buff = webbuff (url, verbose, &lbuff);

	if (buff == NULL) {
	    if (verbose)
		printf ("no return from SIMBAD\n");
	    else {
		if (ned)
		    fprintf (stderr,"*** No return from NED for %s\n",objname);
		if (viz)
		    fprintf (stderr,"*** No return from Vizier for %s\n",objname);
		if (simbad)
		    fprintf (stderr,"*** No return from SIMBAD for %s\n",objname);
		}
	    continue;
	    }
	else if (printall) {
	    buff[lbuff] = (char) 0;
	    xbuff = buff + lbuff - 20;
	    (void) printf ("%s\n", xbuff);
	    (void) printf ("%s\n", buff);
	    (void) printf ("End of Buffer\n");
	    continue;
	    }

	/* Read number of objects identified */
	if ((idline = strsrch (buff, "#=")) != NULL) {
	    id = strchr (idline, ':');
	    if (id != NULL)
		nid = atoi (id+2);
	    else
		nid = 0;

	    /* Get position */
	    if ((posline = strsrch (buff, "%J ")) != NULL) {
		posra = posline + 3;
		while (*posra == ' ')
		    posra++;
		ra = atof (posra);
		posdec = strchr (posra, ' ');
		while (*posdec == ' ')
		    posdec++;
		posline = strchr (posdec, ' ');
		dec = atof (posdec);
		}
	    else {
		if (verbose)
		    printf ("no position from SIMBAD\n");
		else  {
		    if (ned)
		        fprintf (stderr,"*** No NED position for %s\n",objname);
		    if (viz)
		        fprintf (stderr,"*** No Vizier position for %s\n",objname);
		    if (simbad)
		        fprintf (stderr,"*** No SIMBAD position for %s\n",objname);
		    }
		continue;
		}
	    }

	else {
	    nid = 0;
	    if ((errline = strsrch (buff, "#!SIMBAD: ")) != NULL) {
		if ((errend = strchr (errline, '\n')) != NULL)
		   *errend = (char) 0;
		fprintf (stderr, "*** %s\n", errline+10);
		}
	    else {
		fprintf (stderr, "*** No SIMBAD position for %s\n", objname);
		}
	    continue;
	    }

	if (nid > 0) {
	    if (verbose) {
		if (ned && nid == 1)
		    fprintf (stdout, "%d object found in NED: \n", nid);
		else if (ned)
		    fprintf (stdout, "%d objects found in NED: \n", nid);
		else if (viz && nid == 1)
		    fprintf (stdout, "%d object found by VIZIER: \n", nid);
		else if (viz)
		    fprintf (stdout, "%d objects found by VIZIER: \n", nid);
		else if (nid == 1)
		    fprintf (stdout, "%d object found by SIMBAD: \n", nid);
		else
		    fprintf (stdout, "%d objects found by SIMBAD: \n", nid);
		}

	/* Print Starbase header if requested */
	    nobj++;
	    if (tabout && nobj == 1) {
		if (ned)
		    printf ("catalog	NED\n");
		else if (viz)
		    printf ("catalog	VIZIER\n");
		else
		    printf ("catalog	SIMBAD\n");
		if (outsys == WCS_GALACTIC)
		    printf ("radecsys	galactic\n");
		else if (outsys == WCS_ECLIPTIC)
		    printf ("radecsys	ecliptic\n");
		else if (outsys == WCS_B1950)
		    printf ("radecsys	B1950\n");
		else
		    printf ("radecsys	J2000\n");
		if (outsys == WCS_B1950) {
		    printf ("equinox	1950.0\n");
		    printf ("epoch	1950.0\n");
		    }
		else {
		    printf ("equinox	2000.0\n");
		    printf ("epoch	2000.0\n");
		    }
		printf ("program	%s\n",RevMsg);
		if (printname)
		    printf ("id    	");
		printf ("ra      	dec       ");
		printf ("\n");
		if (printname)
		    printf ("------	");
		printf ("------------	------------");
		printf ("\n");
		}
	    if (printname) {
		printf ("%s", objname);
		if (tabout)
		    printf ("	");
		else
		    printf (" ");
		}
	    if (outsys != WCS_J2000)
		wcscon (sysj, outsys, 0.0, 0.0, &ra, &dec, 0.0);
	    if (outsys == WCS_ECLIPTIC || outsys == WCS_GALACTIC) {
		if (verbose)
		    fprintf (stdout, "l= ");
		fprintf (stdout, "%.6f", ra);
		if (tabout)
		    printf ("	");
		else
		    printf (" ");
		if (verbose)
		    fprintf (stdout, "b= ");
		if (dec >= 0.0)
		    fprintf (stdout, "+");
		fprintf (stdout, "%.6f", dec);
		if (!tabout) {
		    if (outsys == WCS_GALACTIC)
			fprintf (stdout, " Galactic");
		    else if (outsys == WCS_ECLIPTIC)
			fprintf (stdout, " Ecliptic");
		    }
		fprintf (stdout, "\n");
		}
	    else {
		if (verbose)
		    fprintf (stdout, "ra= ");
		if (printdeg)
		    fprintf (stdout, "%.6f", ra);
		else {
		    ra2str (rastr, 31, ra, 3);
		    fprintf (stdout, "%s", rastr);
		    }
		if (tabout)
		    printf ("	");
		else
		    printf (" ");
		if (verbose)
		    fprintf (stdout, "dec=");
		if (printdeg)
		    fprintf (stdout, "%.6f", dec);
		else {
		    dec2str (decstr, 31, dec, 2);
		    fprintf (stdout, "%s", decstr);
		    }
		if (!tabout) {
		    if (outsys == WCS_B1950)
			fprintf (stdout, " B1950");
		    else
			fprintf (stdout, " J2000");
		    }
		fprintf (stdout, "\n");
		}
	    if (printid) {
		buffid = buff;
		while ((buffid = strsrch (buffid, "%I ")) != NULL) {
		    buffid = buffid + 3;
		    if ((iend = strchr (buffid, lf))) {
			*iend = (char) 0;
			if ((ieq = strchr (buffid, eq)))
			    *ieq = (char) 0;
			printf ("%s\n", buffid);
			if (ieq)
			    *ieq = eq;
			*iend = lf;
			}
		    else if ((iend = strchr (buffid,cr))) {
			*iend = (char) 0;
			if ((ieq = strchr (buffid, eq)))
			    *ieq = (char) 0;
			printf ("%s\n", buffid);
			if (ieq)
			    *ieq = eq;
			*iend = cr;
			}

		    }
		}
	    }
	free (buff);
	}
    exit (0);
}

static void
PrintUsage (command)

char	*command;	/* Command where error occurred or NULL */

{
    fprintf (stderr,"%s\n", RevMsg);
    fprintf (stderr,"Return RA and Dec for object name using SIMBAD, NED, Vizier\n");
    fprintf (stderr,"Usage:  simpos [-idtv][b|e|g] name1 name2 ...\n");
    fprintf (stderr,"        simpos [-idtv][b|e|g] @namelist ...\n");
    fprintf (stderr,"name(n): Objects for which to search (space -> _)\n");
    fprintf (stderr,"namelist: File with one object name per line\n");
    fprintf (stderr,"-a: Print all information returned\n");
    fprintf (stderr,"-b: Print coordinates in B1950 instead of J2000\n");
    fprintf (stderr,"-c: Use CDS server in France, default is CfA in US\n");
    fprintf (stderr,"-d: Print coordinates in degrees instead of sexigesimal\n");
    fprintf (stderr,"-e: Print coordinates in ecliptic instead of J2000\n");
    fprintf (stderr,"-g: Print coordinates in galactic instead of J2000\n");
    fprintf (stderr,"-i: Print ID(s) returned from SIMBAD\n");
    fprintf (stderr,"-n: Find coordinates using the NED database\n");
    fprintf (stderr,"-s: Find coordinates using SIMBAD\n");
    fprintf (stderr,"-t: Print output as tab-separated table\n");
    fprintf (stderr,"-v: Print extra descriptive info\n");
    fprintf (stderr,"-z: Find coordinates in Vizier catalog\n");
    exit (1);
}

/* Oct 25 2002	New program based on nedpos.c
 *
 * Jun 20 2006	Clean up code
 *
 * Jan 10 2007	Declare RevMsg static, not const
 * Jan 10 2007	exit(0) if successful
 * Jan 11 2007	Include fitsfile.h instead of fitshead.h
 * Jan 16 2007	Fix leading space bug
 * Sep 19 2007	Add -t option for Starbase output and @ for lists of names
 *
 * Apr 16 2008	Fix @ mode for processing lists of objects
 * 
 * Mar  3 2009	Put -i option back in; print name only if requested
 *
 * Jan 28 2013	Update Vizier URL
 * Feb  1 2013	Add -s NSV option to allow users to change search order
 * Feb  1 2013	Add -a option to print everything which is returned
 * Feb  1 2013	Fix -i option to print list of catalog IDs returned
 *
 * May  5 2015	Switch to CfA Vizier server as default
 * May 15 2015	Explicitly search SIMBAD, NED, or Vizier (nedpos or vizpos)
 */

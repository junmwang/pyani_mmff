# include <stdio.h>
# include <math.h>
# include <ctype.h>
# include <stdlib.h>
# include <string.h>
# define MAXCHAR 128
# define COLORTEXT "YES"
# define MAXATOM 512
# define debug 0

typedef struct {
	int  atomic_num;
	char element[5];
	char name[10];
	double x;
	double y;
	double z;
} ATOM;
ATOM atom[MAXATOM];
int natom = 0;
char ifilename[MAXCHAR];
char ofilename[MAXCHAR];
char cfilename[MAXCHAR];
char efilename[MAXCHAR];
char line[MAXCHAR];
int i_input = 0;
int i_output= 0;
int i_cmd = 0;
int i_ene = 0;
int format = 1;
FILE *fpin, *fpout, *fpcmd, *fpene;

int max_force=0;
int max_displacement=0;
int rms_force=0;
int rms_displacement=0;
double energy = 0;

void read() {
int i;
char tmpc1[MAXCHAR];
char tmpc2[MAXCHAR];
char tmpc3[MAXCHAR];
char tmpc4[MAXCHAR];
char tmpc5[MAXCHAR];
char tmpc6[MAXCHAR];
char tmpc7[MAXCHAR];
char tmpc8[MAXCHAR];
char tmpc9[MAXCHAR];
char tmpc10[MAXCHAR];
char tmpc11[MAXCHAR];
for (;;) {
	if (fgets(line, MAXCHAR, fpin) == NULL) {
        	break;
	}
	strcpy(tmpc1, "");
	strcpy(tmpc2, "");
	strcpy(tmpc3, "");
	strcpy(tmpc4, "");
	strcpy(tmpc5, "");
	strcpy(tmpc6, "");
	strcpy(tmpc7, "");
	strcpy(tmpc8, "");
	strcpy(tmpc9, "");
	strcpy(tmpc10, "");
	strcpy(tmpc11, "");

	tmpc1[0]='\0';
	tmpc2[0]='\0';
	tmpc3[0]='\0';
	tmpc4[0]='\0';
	tmpc5[0]='\0';
	tmpc6[0]='\0';
	tmpc7[0]='\0';
	tmpc8[0]='\0';
	tmpc9[0]='\0';
	tmpc10[0]='\0';
	tmpc11[0]='\0';
	sscanf(line, "%s%s%s%s%s%s%s%s%s%s%s", tmpc1, tmpc2, tmpc3, tmpc4, tmpc5, tmpc6, tmpc7, tmpc8, tmpc9, tmpc10, tmpc11);
	if(strcmp(tmpc1, "Maximum") == 0 && strcmp(tmpc2, "Force") == 0) {
		if(strcmp(tmpc5, "YES") == 0) max_force = 1;
	}
	if(strcmp(tmpc1, "RMS") == 0 && strcmp(tmpc2, "Force") == 0) {
		if(strcmp(tmpc5, "YES") == 0) rms_force = 1;
	}
	if(strcmp(tmpc1, "Maximum") == 0 && strcmp(tmpc2, "Displacement") == 0) {
		if(strcmp(tmpc5, "YES") == 0) max_displacement = 1;
	}
	if(strcmp(tmpc1, "RMS") == 0 && strcmp(tmpc2, "Displacement") == 0) {
		if(strcmp(tmpc5, "YES") == 0) rms_displacement = 1;
	}
	if(strcmp(tmpc1, "Energy:") == 0) energy=atof(tmpc2);
	
	if(strcmp(tmpc1, "ATOM") == 0) {
		strcpy(atom[natom].name, tmpc3);
		atom[natom].x= atof(tmpc6);
		atom[natom].y= atof(tmpc7);
		atom[natom].z= atof(tmpc8);
		strcpy(atom[natom].element, tmpc11);
		natom++;
	}
 }
	if(debug) {
		for(i=0; i<natom; i++) 
			printf("%4d %5s %5s %9.4lf %9.4lf %9.4lf\n", i+1, atom[i].name, atom[i].element, atom[i].x, atom[i].y,atom[i].z);
	}
}

void write() {
int i;
char name[10];

for(i=0; i<natom;i++) {
	strcpy(name, atom[i].name);
        if(strlen(name) <= 3)
        	fprintf(fpout, "ATOM%7d  %-3s %3s %5s", i+1, name, "MOL", "1");
        else
                fprintf(fpout, "ATOM%7d %4s %3s %5s", i+1, name, "MOL", "1");
		
	fprintf(fpout, "%12.3lf%8.3lf%8.3lf%6.2lf%6.2lf%12s\n", atom[i].x, atom[i].y, atom[i].z, 1.00, 0.00, atom[i].element);
}
fprintf(fpout, "\n");
if(i_ene == 1) {
	fprintf(fpene, "Energy: %20.12lf\n", energy);
	fprintf(fpene, "MAX_Force converged: %2d\n", max_force);
	fprintf(fpene, "MAX_Displacement: %2d\n", max_displacement);
	fprintf(fpene, "RMS_Force converged: %2d\n", rms_force);
	fprintf(fpene, "RMS_Displacement: %2d\n", rms_displacement);
	if(max_force + max_displacement + rms_force + rms_displacement == 4)
		fprintf(fpene, "Converged: 	Yes\n");
	fclose(fpene);
}
}

void rw_cmd(void) {
int i;
int flag;

char tmpc1[MAXCHAR];
char tmpc2[MAXCHAR];
char tmpc3[MAXCHAR];
char tmpc4[MAXCHAR];
char tmpc5[MAXCHAR];
flag = 0;
for (;;) {
	if (fgets(line, MAXCHAR, fpcmd) == NULL) {
        	break;
	}
	strcpy(tmpc1, "");
	strcpy(tmpc2, "");
	strcpy(tmpc3, "");
	strcpy(tmpc4, "");
	strcpy(tmpc5, "");

	tmpc1[0]='\0';
	tmpc2[0]='\0';
	tmpc3[0]='\0';
	tmpc4[0]='\0';
	tmpc5[0]='\0';
	sscanf(line, "%s%s%s%s%s", tmpc1, tmpc2, tmpc3, tmpc4, tmpc5); 
	if(flag == 0 && strcmp(tmpc1, "#COORD") == 0) {
		flag = 1;
		fprintf(fpout, "%s", line);
		for(i=0;i<natom;i++) 
			fprintf(fpout, "%16.9lf %16.9lf %16.9lf\n", atom[i].x, atom[i].y, atom[i].z);
		fprintf(fpout, "\n");
		continue;
	}
	if(flag == 1 && line[0]=='#')
		flag = 2;
	if(flag == 0 || flag == 2)
		fprintf(fpout, "%s", line);
 }
}

int main(int argc, char *argv[]) {
	int i, flag;
	char cmdstr[MAXCHAR];
	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf
				("[31mUsage: ani2ani -i [0m ani logfile \n"
				 "[31m               -o [0m output file \n"
				 "[31m               -c [0m ani cmd/input file, otional \n"
				 "[31m               -e [0m ene file, optional \n"
				 "[31m               -f [0m output file format, optional \n"
				 "[34m                   1[0m pdb file,the default\n"
				 "[34m                   2[0m new cmd/input file \n" );
			exit(0);
		}
		if (argc != 5 && argc != 7 && argc !=9 && argc != 11 && argc != 13 && argc != 15) {
			printf
				("[31mUsage: ani2ani -i [0m ani logfile \n"
				 "[31m               -o [0m output file \n"
				 "[31m               -c [0m ani cmd/input file, otional \n"
				 "[31m               -e [0m ene file, optional \n"
				 "[31m               -f [0m output file format, optional \n"
				 "[34m                   1[0m pdb file,the default\n"
				 "[34m                   2[0m new cmd/input file \n" );
			exit(0);
		}
	} else {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf
				("Usage: ani2ani -i  ani log file \n"
				 "               -o  output file \n"
				 "               -c  ani cmd/input file, optional\n"
				 "               -e  ene file, optional\n"
				 "               -f  output file format\n"
				 "                   1: pdb file, the default\n" 
				 "                   2: new ani cmd/input\n");
			exit(0);
		}
		if (argc != 5 && argc != 7 && argc !=9 && argc != 11 && argc != 13 && argc != 15) {
			printf
				("Usage: ani2ani -i  ani log file \n"
				 "               -o  output file \n"
				 "               -c  ani cmd/input file, optional\n"
				 "               -e  ene file, optional\n"
				 "               -f  output file format\n"
				 "                   1: pdb file, the default\n" 
				 "                   2: new ani cmd/input\n");
		}
	}

	for (i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-i") == 0) {
			strcpy(ifilename, argv[i + 1]);
			i_input = 1;
		}
		if (strcmp(argv[i], "-o") == 0)  {
			strcpy(ofilename, argv[i + 1]);
			i_output = 1;
		}
		if (strcmp(argv[i], "-c") == 0) {
			strcpy(cfilename, argv[i + 1]);
			i_cmd = 1;
		}
		if (strcmp(argv[i], "-e") == 0) {
			strcpy(efilename, argv[i + 1]);
			i_ene = 1;
		}
		if (strcmp(argv[i], "-f") == 0) {
			format=atoi(argv[i + 1]);
		}
	}
	if (format != 1 && format != 2) format = 1;
	if (format == 2 && i_cmd == 0) {
		fprintf(stderr, "To output a new cmd/input file, a cmd/input file must be provided with '-c' flag\n");
		exit(1);
	}
	if ((fpin = fopen(ifilename, "r")) == NULL) {
        	fprintf(stderr, "Cannot open ani input file %s to read, exit\n", ifilename);
        	exit(1);
	}
	if ((fpout = fopen(ofilename, "w")) == NULL) {
        	fprintf(stderr, "Cannot open output file %s to write, exit\n", ofilename);
        	exit(1);
	}
	if(i_ene==1) { 
		if ((fpene = fopen(efilename, "w")) == NULL) {
        		fprintf(stderr, "Cannot open energy file %s to write, exit\n", efilename);
        		exit(1);
		}
	}
	if(i_cmd==1 && format ==2) { 
		if ((fpcmd = fopen(cfilename, "r")) == NULL) {
        		fprintf(stderr, "Cannot open ani log file %s to read, exit\n", cfilename);
        		exit(1);
		}
	}
	read();
	if(format == 1) write();
	if(i_cmd  == 1 && format == 2) 
		rw_cmd();
fclose(fpin);
fclose(fpout);
return 0;
}

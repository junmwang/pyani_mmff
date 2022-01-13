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
	double x;
	double y;
	double z;
} ATOM;
ATOM atom[MAXATOM];
int natom = 0;
char ifilename[MAXCHAR];
char ofilename[MAXCHAR];
char line[MAXCHAR];
FILE *fpin, *fpout; 
int i_input = 0;
int i_output= 0;
int i_cb    = 0;
int i_ca    = 0;
int i_ct    = 0;
int i_opt   = 0;
int i_comment = 0;
int maxstep = -1;
char cb_str[MAXCHAR]="";
char ca_str[MAXCHAR]="";
char ct_str[MAXCHAR]="";
char compname[MAXCHAR]="MOL";
char comment[MAXCHAR];
char elem[MAXCHAR];
int charge = 0;
int multiplicity = 1;

void element(int atomicnum) {
	int i;
	switch (atomicnum) {
	case 1:
		strcpy(elem, "H");
		break;
	case 2:
		strcpy(elem, "He");
		break;
	case 3:
		strcpy(elem, "Li");
		break;
	case 4:
		strcpy(elem, "Be");
		break;
	case 5:
		strcpy(elem, "B");
		break;
	case 6:
		strcpy(elem, "C");
		break;
	case 7:
		strcpy(elem, "N");
		break;
	case 8:
		strcpy(elem, "O");
		break;
	case 9:
		strcpy(elem, "F");
		break;
	case 10:
		strcpy(elem, "Ne");
		break;
	case 11:
		strcpy(elem, "Na");
		break;
	case 12:
		strcpy(elem, "Mg");
		break;
	case 13:
		strcpy(elem, "Al");
		break;
	case 14:
		strcpy(elem, "Si");
		break;
	case 15:
		strcpy(elem, "P");
		break;
	case 16:
		strcpy(elem, "S");
		break;
	case 17:
		strcpy(elem, "Cl");
		break;
	case 18:
		strcpy(elem, "Ar");
		break;
	case 19:
		strcpy(elem, "K");
		break;
	case 20:
		strcpy(elem, "Ca");
		break;
	case 21:
		strcpy(elem, "Sc");
		break;
	case 22:
		strcpy(elem, "Ti");
		break;
	case 23:
		strcpy(elem, "V");
		break;
	case 24:
		strcpy(elem, "Cr");
		break;
	case 25:
		strcpy(elem, "Mn");
		break;
	case 26:
		strcpy(elem, "Fe");
		break;
	case 27:
		strcpy(elem, "Co");
		break;
	case 28:
		strcpy(elem, "Ni");
		break;
	case 29:
		strcpy(elem, "Cu");
		break;
	case 30:
		strcpy(elem, "Zn");
		break;
	case 31:
		strcpy(elem, "Ga");
		break;
	case 32:
		strcpy(elem, "Ge");
		break;
	case 33:
		strcpy(elem, "As");
		break;
	case 34:
		strcpy(elem, "Se");
		break;
	case 35:
		strcpy(elem, "Br");
		break;
	case 36:
		strcpy(elem, "Kr");
		break;
	case 37:
		strcpy(elem, "Rb");
		break;
	case 38:
		strcpy(elem, "Sr");
		break;
	case 39:
		strcpy(elem, "Y");
		break;
	case 40:
		strcpy(elem, "Zr");
		break;
	case 41:
		strcpy(elem, "Nb");
		break;
	case 42:
		strcpy(elem, "Mo");
		break;
	case 43:
		strcpy(elem, "Tc");
		break;
	case 44:
		strcpy(elem, "Ru");
		break;
	case 45:
		strcpy(elem, "Rh");
		break;
	case 46:
		strcpy(elem, "Pd");
		break;
	case 47:
		strcpy(elem, "Ag");
		break;
	case 48:
		strcpy(elem, "Cd");
		break;
	case 49:
		strcpy(elem, "In");
		break;
	case 50:
		strcpy(elem, "Sn");
		break;
	case 51:
		strcpy(elem, "Sb");
		break;
	case 52:
		strcpy(elem, "Te");
		break;
	case 53:
		strcpy(elem, "I");
		break;
	case 54:
		strcpy(elem, "Xe");
		break;
	case 55:
		strcpy(elem, "Cs");
		break;
	case 56:
		strcpy(elem, "Ba");
		break;
	case 57:
		strcpy(elem, "La");
		break;
	case 58:
		strcpy(elem, "Ce");
		break;
	case 59:
		strcpy(elem, "Pr");
		break;
	case 60:
		strcpy(elem, "Nd");
		break;
	case 61:
		strcpy(elem, "Pm");
		break;
	case 62:
		strcpy(elem, "Sm");
		break;
	case 63:
		strcpy(elem, "Eu");
		break;
	case 64:
		strcpy(elem, "Gd");
		break;
	case 65:
		strcpy(elem, "Tb");
		break;
	case 66:
		strcpy(elem, "Dy");
		break;
	case 67:
		strcpy(elem, "Ho");
		break;
	case 68:
		strcpy(elem, "Er");
		break;
	case 69:
		strcpy(elem, "Tm");
		break;
	case 70:
		strcpy(elem, "Yb");
		break;
	case 71:
		strcpy(elem, "Lu");
		break;
	case 72:
		strcpy(elem, "Hf");
		break;
	case 73:
		strcpy(elem, "Ta");
		break;
	case 74:
		strcpy(elem, "W");
		break;
	case 75:
		strcpy(elem, "Re");
		break;
	case 76:
		strcpy(elem, "Os");
		break;
	case 77:
		strcpy(elem, "Ir");
		break;
	case 78:
		strcpy(elem, "Pt");
		break;
	case 79:
		strcpy(elem, "Au");
		break;
	case 80:
		strcpy(elem, "Hg");
		break;
	case 81:
		strcpy(elem, "Tl");
		break;
	case 82:
		strcpy(elem, "Pb");
		break;
	case 83:
		strcpy(elem, "Bi");
		break;
	case 84:
		strcpy(elem, "Po");
		break;
	case 85:
		strcpy(elem, "At");
		break;
	case 86:
		strcpy(elem, "Rn");
		break;
	case 87:
		strcpy(elem, "Fr");
		break;
	case 88:
		strcpy(elem, "Ra");
		break;
	case 89:
		strcpy(elem, "Ac");
		break;
	case 90:
		strcpy(elem, "Th");
		break;
	case 91:
		strcpy(elem, "Pa");
		break;
	case 92:
		strcpy(elem, "U");
		break;
	case 93:
		strcpy(elem, "Np");
		break;
	case 94:
		strcpy(elem, "Pu");
		break;
	case 95:
		strcpy(elem, "Am");
		break;
	case 96:
		strcpy(elem, "Cm");
		break;
	case 97:
		strcpy(elem, "Bk");
		break;
	case 98:
		strcpy(elem, "Cf");
		break;
	case 99:
		strcpy(elem, "Es");
		break;
	case 100:
		strcpy(elem, "Fm");
		break;
	case 101:
		strcpy(elem, "Md");
		break;
	case 102:
		strcpy(elem, "No");
		break;
	case 103:
		strcpy(elem, "Lr");
		break;
	case 104:
		strcpy(elem, "Rf");
		break;
	case 105:
		strcpy(elem, "Db");
		break;
	case 106:
		strcpy(elem, "Sg");
		break;
	case 107:
		strcpy(elem, "Bh");
		break;
	case 108:
		strcpy(elem, "Hs");
		break;
	case 109:
		strcpy(elem, "Mt");
		break;
	case 110:
		strcpy(elem, "Ds");
		break;
	default:
		strcpy(elem, "du");
		break;
	}
}

int atomic_num(char *elem) {
int atomicnum;
	switch (elem[0]) {
		case 'A':
			if (elem[1] == 'c') 
				atomicnum = 89;
			else if (elem[1] == 'g') 
				atomicnum = 47;
			else if (elem[1] == 'l')
				atomicnum = 13;
			else if (elem[1] == 'm')
				atomicnum = 95;
			else if (elem[1] == 'r') 
				atomicnum = 18;
			else if (elem[1] == 's') 
				atomicnum = 33;
			else if (elem[1] == 't') 
				atomicnum = 85;
			else if (elem[1] == 'u') 
				atomicnum = 79;
			break;
		case 'B':
			if (elem[1] == 'r' || elem[1] == 'R') 
				atomicnum = 35;
			else if (elem[1] == 'a') 
				atomicnum = 56;
			else if (elem[1] == 'e') 
				atomicnum = 4;
			else if (elem[1] == 'h') 
				atomicnum = 107;
			else if (elem[1] == 'i') 
				atomicnum = 83;
			else if (elem[1] == 'k') 
				atomicnum = 97;
			else 
				atomicnum = 5;
			break;
		case 'C':
			if (elem[1] == 'l' || elem[1] == 'L') 
				atomicnum = 17;
			else if (elem[1] == 'a') 
				atomicnum = 20;
			else if (elem[1] == 'd')  
				atomicnum = 48;
			else if (elem[1] == 'e')  
				atomicnum = 58;
			else if (elem[1] == 'f')  
				atomicnum = 98;
			else if (elem[1] == 'm')  
				atomicnum = 96;
			else if (elem[1] == 'o') 
				atomicnum = 27;
			else if (elem[1] == 'r')
				atomicnum = 24;
			else if (elem[1] == 's')
				atomicnum = 55;
			else if (elem[1] == 'u')
				atomicnum = 29;
			else 
				atomicnum = 6;
			break;
		case 'D':
			if (elem[1] == 'b')  
				atomicnum = 105;
			else if (elem[1] == 's')
				atomicnum = 110;
			else if (elem[1] == 'y') 
				atomicnum = 66;
			else 
				atomicnum = 1;
			break;
		case 'E':
			if (elem[1] == 'P') 
				atomicnum = 0;
			else if (elem[1] == 'r')
				atomicnum = 68;
			else if (elem[1] == 's')
				atomicnum = 99;
			else if (elem[1] == 'u')
				atomicnum = 63;
			break;
		case 'F':
			if (elem[1] == 'e') 
				atomicnum = 26;
			else if (elem[1] == 'm')
				atomicnum = 100;
			else if (elem[1] == 'r')
				atomicnum = 87;
			else
				atomicnum = 9;
			break;
		case 'G':
			if (elem[1] == 'a') 
				atomicnum = 31;
			else if (elem[1] == 'd')
				atomicnum = 64;
			else if (elem[1] == 'e')
				atomicnum = 32;
			break;
		case 'H':
			if (elem[1] == 'e')
				atomicnum = 2;
			else if (elem[1] == 'f')
				atomicnum = 72;
			else if (elem[1] == 'g') 
				atomicnum = 80;
			else if (elem[1] == 'o')
				atomicnum = 67;
			else if (elem[1] == 's') 
				atomicnum = 108;
			else 
				atomicnum = 1;
			break;
		case 'I':
			if (elem[1] == 'n')
				atomicnum = 49;
			else if (elem[1] == 'r')
				atomicnum = 77;
			else 
				atomicnum = 53;
			break;
		case 'K':
			if (elem[1] == 'r') 
				atomicnum = 36;
			else 
				atomicnum = 19;
			break;
		case 'l':
			if (elem[1] == 'p')
				atomicnum = 0;
			break;
		case 'L':
			if (elem[1] == 'i') 
				atomicnum = 3;
			else if (elem[1] == 'a') 
				atomicnum = 57;
			else if (elem[1] == 'r') 
				atomicnum = 103;
			else if (elem[1] == 'u') 
				atomicnum = 71;
			else if (elem[1] == 'P') 
				atomicnum = 0;
			break;
		case 'M':
			if (elem[1] == 'n') 
				atomicnum = 25;
			else if (elem[1] == 'g') 
				atomicnum = 12;
			else if (elem[1] == 'd')
				atomicnum = 101;
			else if (elem[1] == 'o')
				atomicnum = 42;
			else if (elem[1] == 't') 
				atomicnum = 109;
			break;
		case 'N':
			if (elem[1] == 'i')
				atomicnum = 28;
			else if (elem[1] == 'a') 
				atomicnum = 11;
			else if (elem[1] == 'b')
				atomicnum = 41;
			else if (elem[1] == 'd')
				atomicnum = 60;
			else if (elem[1] == 'e')
				atomicnum = 10;
			else if (elem[1] == 'o') 
				atomicnum = 102;
			else if (elem[1] == 'p')
				atomicnum = 93;
			else 
				atomicnum = 7;
			break;
		case 'O':
			if (elem[1] == 's') 
				atomicnum = 76;
			else
				atomicnum = 8;
			break;
		case 'P':
			if (elem[1] == 'd') 
				atomicnum = 46;
			else if (elem[1] == 't') 
				atomicnum = 78;
			else if (elem[1] == 'b') 
				atomicnum = 82;
			else if (elem[1] == 'a')
				atomicnum = 91;
			else if (elem[1] == 'm')
				atomicnum = 61;
			else if (elem[1] == 'o')
				atomicnum = 84;
			else if (elem[1] == 'r')
				atomicnum = 59;
			else if (elem[1] == 'u') 
				atomicnum = 94;
			else 
				atomicnum = 15;
			break;
		case 'R':
			if (elem[1] == 'u')
				atomicnum = 44;
			else if (elem[1] == 'h')
				atomicnum = 45;
			else if (elem[1] == 'a')
				atomicnum = 88;
			else if (elem[1] == 'b')
				atomicnum = 37;
			else if (elem[1] == 'e')
				atomicnum = 75;
			else if (elem[1] == 'f')
				atomicnum = 104;
			else if (elem[1] == 'n') 
				atomicnum = 86;
			break;
		case 'S':
			if (elem[1] == 'i' || elem[1] == 'I')
				atomicnum = 14;
			else if (elem[1] == 'c')
				atomicnum = 21;
			else if (elem[1] == 'e') 
				atomicnum = 34;
			else if (elem[1] == 'r') 
				atomicnum = 38;
			else if (elem[1] == 'b') 
				atomicnum = 51;
			else if (elem[1] == 'g')
				atomicnum = 106;
			else if (elem[1] == 'm')
				atomicnum = 62;
			else if (elem[1] == 'n')
				atomicnum = 50;
			else  
				atomicnum = 16;
			break;
		case 'T':
			if (elem[1] == 'i') 
				atomicnum = 22;
			else if (elem[1] == 'l') 
				atomicnum = 81;
			else if (elem[1] == 'a')
				atomicnum = 73;
			else if (elem[1] == 'b') 
				atomicnum = 65;
			else if (elem[1] == 'c') 
				atomicnum = 43;
			else if (elem[1] == 'e') 
				atomicnum = 52;
			else if (elem[1] == 'h')
				atomicnum = 90;
			else if (elem[1] == 'm')
				atomicnum = 69;
			else 
				atomicnum = 1;
			break;
		case 'U':
			atomicnum = 92;
			break;
		case 'V':
			atomicnum = 23;
			break;
		case 'W':
			atomicnum = 74;
			break;
		case 'X':
			if (elem[1] == 'e') 
				atomicnum = 54;
			break;
		case 'Y':
			if (elem[1] == 'b') 
				atomicnum = 70;
			else 
				atomicnum = 39;
			break;
		case 'Z':
			if (elem[1] == 'n') 
				atomicnum = 30;
			else if (elem[1] == 'r') 
				atomicnum = 40;
			break;
		default:
			printf("\n Unrecognized atomic elem %5s, exit", elem);
	}
	return atomicnum;
}


void read() {
int i;
int read_flag = 0;
int nspaceline=0;
int len;
int optflag = 0;
int tmpint;
char tmpc1[MAXCHAR];
char tmpc2[MAXCHAR];
char tmpc3[MAXCHAR];
char tmpc4[MAXCHAR];
char tmpc5[MAXCHAR];

for (;;) {
	if (fgets(line, MAXCHAR, fpin) == NULL) {
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
	len=strlen(tmpc1);
	if(len == 0) nspaceline++;
        if(nspaceline == 0 && read_flag == 0) {
		if(strncmp(line, "%chk=", 5) == 0) {
			for(i=5; i<=len; i++)
				compname[i-5]=line[i];	
		}
                if(line[0]=='#') {
                        if(strstr(line, "opt") != NULL)
                                i_opt=1;
                        if(strstr(line, "Opt") != NULL)
                                i_opt=1;
                        if(strstr(line, "OPT") != NULL)
                                i_opt=1;
                }
        }
	if(nspaceline == 1 && read_flag == 0 && i_comment == 0) {
		fgets(line, MAXCHAR, fpin); 
		i_comment = 1;
		strcpy(comment, line);
	}
	if(nspaceline == 2 && read_flag == 0) {
		fgets(line, MAXCHAR, fpin); 
		sscanf(line, "%d%d", &charge, &multiplicity);
		read_flag = 1;
		continue;
	}
	if(read_flag == 1 && len == 0) {
		read_flag = -1;
		break;
	}
	if(read_flag == 1) {
		if(optflag == 0 && strlen(tmpc5) > 1) optflag = 1;
		if(isdigit(tmpc1[0]) != 0) {
			if(optflag == 0) sscanf(line, "%d %lf %lf %lf\n", &atom[natom].atomic_num, &atom[natom].x, &atom[natom].y, &atom[natom].z);
			if(optflag == 1) sscanf(line, "%d %d %lf %lf %lf\n", &atom[natom].atomic_num, &tmpint, &atom[natom].x, &atom[natom].y, &atom[natom].z);
			strcpy(atom[natom].element, "");
		}
		else {
			if(optflag == 0) sscanf(line, "%s %lf %lf %lf\n", atom[natom].element, &atom[natom].x, &atom[natom].y, &atom[natom].z);
			if(optflag == 1) sscanf(line, "%s %d %lf %lf %lf\n", atom[natom].element, &tmpint, &atom[natom].x, &atom[natom].y, &atom[natom].z);
			atom[natom].atomic_num = -1;
		}
			natom++;
	}
   }
	for(i=0; i<natom ; i++) {
		if(atom[i].atomic_num == -1) atom[i].atomic_num = atomic_num(atom[i].element);
		if(atom[i].atomic_num != -1) {
				strcpy(elem, "du");
				element(atom[i].atomic_num);
				strcpy(atom[i].element, elem);
		}
	}
	if(debug) {
		for(i=0; i<natom; i++) 
			printf("%4d %5d %5s %9.4lf %9.4lf %9.4lf\n", i+1, atom[i].atomic_num, atom[i].element, atom[i].x, atom[i].y,atom[i].z);
	}
}

void write() {
int i;
int count;
char str[MAXCHAR];

fprintf(fpout, "#MOLNAME\n%s\n", compname);
if(i_opt==1)
	fprintf(fpout, "#MAXSTEP\n%d\n", maxstep);
else
	fprintf(fpout, "#MAXSTEP\n%d\n",  1);
fprintf(fpout, "! -1 for using the default value\n\n");
fprintf(fpout, "#CONVERGE\n%d\n", 1);
fprintf(fpout, "! 1: using the gaussian default criteria\n");
fprintf(fpout, "! 2:0.5, 3:0.6, 4:0.7, 5:0.8, 6:0.9 -using more stringent criteria\n");
fprintf(fpout, "! 7:1.1, 4:1.2, 5:1.3, 6:1.4, 7:1.5, 8:1.6, 9:1.7, 10:1.8 -using less stringent criteria\n");
fprintf(fpout, "! 11:2.0, 12:2.2, 13:2.4, 14:2.6, 15:2.8, 16:3.0, 17:4.0, 18:5.0 -using losse criteria\n");
fprintf(fpout, "\n");
fprintf(fpout, "#COMMENT\n%s\n", comment);
fprintf(fpout, "#CHARGE\n%d\n\n", charge);
fprintf(fpout, "#MULTI\n%d\n\n", multiplicity);
fprintf(fpout, "#ELEM\n");
count = 0;
for(i=0; i<natom;i++) {
	fprintf(fpout, "%3s ", atom[i].element);
	count++;
	if(count == 10) {
		count = 0;
		fprintf(fpout, "\n");		
	}
}
if(count < 10) 
	fprintf(fpout, "\n\n");
else
	fprintf(fpout, "\n");
fprintf(fpout, "#ATOMIC_NUM\n");
count = 0;
for(i=0; i<natom;i++) {
	fprintf(fpout, "%3d ", atom[i].atomic_num);
	count++;
	if(count == 10) {
		count = 0;
		fprintf(fpout, "\n");		
	}
}
if(count < 10) 
	fprintf(fpout, "\n\n");
else
	fprintf(fpout, "\n");
fprintf(fpout, "#COORD\n");
for(i=0; i<natom;i++) {
	fprintf(fpout, "%16.9lf %16.9lf %16.9lf\n", atom[i].x, atom[i].y, atom[i].z);
}
fprintf(fpout, "\n");
fprintf(fpout, "#CONS_BOND\n");
if(i_cb == 1) {
	count = 0;
	for(i=0; i<strlen(cb_str); i++) {
		if(cb_str[i]=='"') {
			count++;
			if(count == 1) continue;		
			if(count == 2) break;
		}
		if(cb_str[i] == '-') cb_str[i] = ' ';
		if(cb_str[i] == ',') {
			 fprintf(fpout, "\n");
			 continue;
		}
		fprintf(fpout, "%c", cb_str[i]);
	}	
	fprintf(fpout, "\n");
}
fprintf(fpout, "\n");


fprintf(fpout, "#CONS_ANGLE\n");
if(i_ca == 1) {
        count = 0;
        for(i=0; i<strlen(ca_str); i++) {
                if(ca_str[i]=='"') {
                        count++;
                        if(count == 1) continue;
                        if(count == 2) break;
                }
                if(ca_str[i] == '-') ca_str[i] = ' ';
                if(ca_str[i] == ',') {
                         fprintf(fpout, "\n");
                         continue;
                }
                fprintf(fpout, "%c", ca_str[i]);
        }
        fprintf(fpout, "\n");
}
fprintf(fpout, "\n");

fprintf(fpout, "#CONS_TOR\n");
if(i_ct == 1) {
        count = 0;
        for(i=0; i<strlen(ct_str); i++) {
                if(ct_str[i]=='"') {
                        count++;
                        if(count == 1) continue;
                        if(count == 2) break;
                }
                if(ct_str[i] == '-') ct_str[i] = ' ';
                if(ct_str[i] == ',') {
                         fprintf(fpout, "\n");
                         continue;
                }
                fprintf(fpout, "%c", ct_str[i]);
        }
        fprintf(fpout, "\n");
}
fprintf(fpout, "\n");


fprintf(fpout, "#END\n");
}

int main(int argc, char *argv[]) {
	int i, flag;
	char cmdstr[MAXCHAR];
	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf
				("[31mUsage: gcrt2ani -i [0m input file \n"
				 "[31m                -o [0m output file \n"
				 "[31m                -ns[0m number of steps\n"
				 "[34m                   [0m -1: using the default nstep value-10000, the default\n" 
				 "[34m                   [0m  n: the maximum step is n\n" 
				 "[31m                -cc[0m convergen creteria, -1 (using the default values)), must be an positive integer\n"
				 "[34m                   [0m -1: using the default gaussian convergence creteria, the default\n" 
				 "[34m                   [0m  v: multiple default values by v\n" 
				 "[31m                -cb[0m constraint string for bond, inside \"\", format: at1-at2, seperated by ',', optional\n" 
				 "[31m                -ca[0m constraint string for angle, inside \"\", format: at1-at2, seperated by ',', optional\n" 
				 "[31m                -ct[0m constraint string for torsion, inside \"\", format: at1-at2, seperated by ',', optional\n" );
			exit(0);
		}
		if (argc != 5 && argc != 7 && argc !=9 && argc != 11 && argc != 13 && argc != 15) {
			printf
				("[31mUsage: gcrt2ani -i [0m input file \n"
				 "[31m                -o [0m output file \n"
				 "[31m                -ns[0m number of steps\n"
				 "[34m                   [0m -1: using the default nstep value-10000, the default\n" 
				 "[34m                   [0m  n: the maximum step is n\n" 
				 "[31m                -cc[0m convergen creteria, -1 (using the default values)), must be an positive integer\n"
				 "[34m                   [0m -1: using the default gaussian convergence creteria, the default\n" 
				 "[34m                   [0m  v: multiple default values by v\n" 
				 "[31m                -cb[0m constraint string for bond, inside \"\", format: at1-at2, seperated by ',', optional\n" 
				 "[31m                -ca[0m constraint string for angle, inside \"\", format: at1-at2, seperated by ',', optional\n" 
				 "[31m                -ct[0m constraint string for torsion, inside \"\", format: at1-at2, seperated by ',', optional\n" );
			exit(0);
		}
	} else {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf
				("Usage: gcrt2ani -i  input file \n"
				 "                -o  output file \n"
				 "                -ns number of steps\n"
				 "                    -1: using the default nstep value-10000, the default\n" 
				 "                     n: the maximum step is n\n" 
				 "                -cc convergen creteria, -1 (using the default values)), must be an positive integer\n"
				 "                    -1: using the default gaussian convergence creteria, the default\n" 
				 "                     v: multiple default values by v\n" 
				 "                -cb constraint string for bond, inside \"\", format: at1-at2, seperated by ',', optional\n" 
				 "                -ca constraint string for angle, inside \"\", format: at1-at2, seperated by ',', optional\n" 
				 "                -ct constraint string for torsion, inside \"\", format: at1-at2, seperated by ',', optional\n" );
			exit(0);
		}
		if (argc != 5 && argc != 7 && argc !=9 && argc != 11 && argc != 13 && argc != 15) {
			printf
				("Usage: gcrt2ani -i  input file \n"
				 "                -o  output file \n"
				 "                -cb constraint string for bond, inside \"\", format: at1-at2, seperated by ',', optional\n" 
				 "                -ca constraint string for angle, inside \"\", format: at1-at2, seperated by ',', optional\n" 
				 "                -ct constraint string for torsion, inside \"\", format: at1-at2, seperated by ',', optional\n" );
			exit(0);
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
		if (strcmp(argv[i], "-cb") == 0) {
			strcpy(cb_str, argv[i + 1]);
			i_cb = 1;
		}
		if (strcmp(argv[i], "-ca") == 0) {
			strcpy(ca_str, argv[i + 1]);
			i_ca = 1;
		}
		if (strcmp(argv[i], "-ct") == 0) {
			strcpy(ct_str, argv[i + 1]);
			i_ct = 1;
		}
		if (strcmp(argv[i], "-ns") == 0) {
			maxstep=atoi(argv[i+1]);
		}
	}
	if ((fpin = fopen(ifilename, "r")) == NULL) {
        	fprintf(stderr, "Cannot open pdb file %s to read, exit\n", ifilename);
        	exit(1);
	}
	if ((fpout = fopen(ofilename, "w")) == NULL) {
        	fprintf(stderr, "Cannot open cmd file %s to write, exit\n", ofilename);
        	exit(1);
	}
	read();
	write();
/*
        fprintf(fpout, "ATOM%7d %-4s %3s  %4d    %8.3f%8.3f%8.3f%6.2lf%6.2lf%12s\n", i + 1, atom_name, "MOL", 1, X, Y, Z, 1.0, 0.0, element); 
*/
fclose(fpin);
fclose(fpout);
return 0;
}

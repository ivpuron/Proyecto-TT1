#include "..\include\globals.h"


Matrix *Global::eopdata;
Matrix *Global::Cnm;
Matrix *Global::Snm;
Matrix *Global::PC;
Matrix *Global::obs;
Param Global::auxparam;


void Global::eop19620101(int c){
    Global::eopdata = new Matrix(13,c);
    FILE *fid = fopen("../data/eop19620101.txt","r");
    if(fid == NULL){
        printf("ERROR");
        exit(EXIT_FAILURE);
    }

    for(int i=1;i<=c;i++){
        fscanf(fid,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
                                                                            &((*Global::eopdata)(1,i)),
                                                                            &((*Global::eopdata)(2,i)),
                                                                            &((*Global::eopdata)(3,i)),
                                                                            &((*Global::eopdata)(4,i)),
                                                                            &((*Global::eopdata)(5,i)),
                                                                            &((*Global::eopdata)(6,i)),
                                                                            &((*Global::eopdata)(7,i)),
                                                                            &((*Global::eopdata)(8,i)),
                                                                            &((*Global::eopdata)(9,i)),
                                                                            &((*Global::eopdata)(10,i)),
                                                                            &((*Global::eopdata)(11,i)),
                                                                            &((*Global::eopdata)(12,i)),
                                                                            &((*Global::eopdata)(13,i)));
        fclose(fid);
    }
}

void Global::GGM03S()
{
    Global::Cnm = new Matrix(181, 181);
    Global::Snm = new Matrix(181, 181);

    FILE *fid = fopen("../data/GGM03S.txt", "r");
    if (fid == NULL)
    {
        printf("Error al abrir el archivo");
        exit(EXIT_FAILURE);
    }

    Matrix aux(6);

    for (int n = 0; n <= 180; n++)
    {
        for (int m = 0; m <= n; m++)
        {
            fscanf(fid, "%lf %lf %lf %lf %lf %lf",
                   &(aux(1)),
                   &(aux(2)),
                   &(aux(3)),
                   &(aux(4)),
                   &(aux(5)),
                   &(aux(6)));

            (*Global::Cnm)(n + 1, m + 1) = aux(3);
            (*Global::Snm)(n + 1, m + 1) = aux(4);
        }
    }

    fclose(fid);
}

void Global::DE430Coeff(int f, int c)
{
    Global::PC = new Matrix(f, c);

    FILE *fid = fopen("../data/DE430Coeff.txt", "r");
    if (fid == NULL)
    {
        printf("Error al abrir el archivo");
        exit(EXIT_FAILURE);
    }

    for (int i = 1; i <= f; i++)
    {
        for (int j = 1; j <= c; j++)
        {
            fscanf(fid, "%lf", &((*Global::PC)(i, j)));
        }
    }

    fclose(fid);
}

void Global::GEOS3(int f){
    Global::obs = new Matrix(f, 4);

    FILE *fid=fopen("../Data/GEOS3.txt","r");
	if(fid==NULL){
		printf("Error al abrir GEOS3.txt\n");
		exit(EXIT_FAILURE);
	}

    int Y, MO, D, H, MI, S;
    double AZ, EL, DIST;

    char line[55], y[5], mo[3], d[3], h[3], mi[3], s[7], az[9], el[9], dist[10];

    for(int i=1; i<=f; i++){
        fgets(line,sizeof(line)+2, fid);

        strncpy(y, &line[0], 4);
        y[4] = '\0';

        strncpy(mo, &line[5], 2);
        mo[2] = '\0';

        strncpy(d, &line[8], 2);
        d[2] = '\0';

        strncpy(h, &line[12], 2);
        h[2] = '\0';

        strncpy(mi, &line[15], 2);
        mi[2] = '\0';

        strncpy(s, &line[18], 6);
        s[6] = '\0';

        strncpy(az, &line[25], 8);
        az[8] = '\0';

        strncpy(el, &line[35], 8);
        el[8] = '\0';

        strncpy(dist, &line[44], 9);
        dist[9] = '\0';

        Y = atoi(y);
        MO = atoi(mo);
        D = atoi(d);
        H = atoi(h);
        MI = atoi(mi);
        S = atoi(s);
        AZ = atof(az);
        EL = atof(el);
        DIST = atof(dist);

        (*Global::obs)(i, 1) = Mjday(Y,MO,D,H,MI,S);
		(*Global::obs)(i, 2) = Const::Rad*AZ;
		(*Global::obs)(i, 3) = Const::Rad*EL;
		(*Global::obs)(i, 4) = 1e3*DIST;

    }
    fclose(fid);
}
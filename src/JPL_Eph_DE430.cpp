


#include "../include/JPL_Eph_DE430.h"

void JPL_Eph_DE430(Matrix & r_Mercury,Matrix & r_Venus,Matrix & r_Earth,Matrix & r_Mars,Matrix & r_Jupiter,Matrix & r_Saturn,Matrix & r_Uranus, 
          Matrix & r_Neptune,Matrix & r_Pluto,Matrix & r_Moon, Matrix & r_Sun, double Mjd_TDB){

    double JD = Mjd_TDB + 2400000.5;
    int i;
    for (int m = 1; m <= 2285; m++)
    {
        if ((*Global::PC)(m, 1) <= JD && JD <= (*Global::PC)(m, 2))
        {
            i=m;
            break;
        }
    }
    Matrix PCtemp = (*Global::PC).columna(i);

    double t1 = PCtemp(1)-2400000.5; // MJD at start of interval

    double dt = Mjd_TDB - t1;

    //temp = (231:13:270);
    Matrix temp(4);
    int ii=1;
    for(int k=231; k<=270; k=k+13){
        temp(ii)=k;
        i++;
    }
    Matrix Cx_Earth (26);
    Matrix Cy_Earth (26);
    Matrix Cz_Earth (26);

    for (int n = 1; n <= 13; n++)
    {
        Cx_Earth(n) = PCtemp(temp(1) + n - 1);
        Cy_Earth(n) = PCtemp(temp(2) + n - 1);
        Cz_Earth(n) = PCtemp(temp(3) + n - 1);
    }

    temp = temp+39;

    Matrix Cx (13);
    Matrix Cy (13);
    Matrix Cz (13);

    for (int n = 1; n <= 13; n++)
    {
        Cx(n) = PCtemp(temp(1) + n - 1);
        Cy(n) = PCtemp(temp(2) + n - 1);
        Cz(n) = PCtemp(temp(3) + n - 1);
    } 

    for (int n = 1; n <= 13; n++)
    {
        Cx_Earth(n + 13) = Cx(n);
        Cy_Earth(n + 13) = Cy(n);
        Cz_Earth(n + 13) = Cz(n);
    }

    int j=-1;
    double Mjd0=0.0;
    if (0<=dt && dt<=16){
        j=0;
        Mjd0 = t1;
    }
    else if(16<dt && dt<=32){
        j=1;
        Mjd0 = t1+16*j;
    }

    Matrix cx(13);
    Matrix cy(13);
    Matrix cz(13);

    for (int n = 1; n <= 13; n++)
    {
        cx(n) = Cx_Earth(13 * j + n);
        cy(n) = Cy_Earth(13 * j + n);
        cz(n) = Cz_Earth(13 * j + n);
    }

    r_Earth = transpose(Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+16, cx,
                        cy, cz))*1e3;

    //temp = (441:13:480);

     ii = 1;
    for (int k = 441; k <= 480; k = k + 13)
    {
        temp(ii) = k;
        i++;
    }

     Matrix Cx_Moon(109);
    Matrix Cy_Moon(109);
    Matrix Cz_Moon(109);
    for (int n = 1; n <= 13; n++)
    {
        Cx_Moon(n) = PCtemp(temp(1) + n - 1);
        Cy_Moon(n) = PCtemp(temp(2) + n - 1);
        Cz_Moon(n) = PCtemp(temp(3) + n - 1);
    }
    for(int i=1;i<=7;i++){
        temp = temp+39;
        for (int n = 1; n <= 13; n++)
        {
            Cx(n) = PCtemp(temp(1) + n - 1);
            Cy(n) = PCtemp(temp(2) + n - 1);
            Cz(n) = PCtemp(temp(3) + n - 1);
        }

        for (int n = 1; n <= 13; n++)
        {
            Cx_Moon(n + 13*i) = Cx(n);
            Cy_Moon(n + 13*i) = Cy(n);
            Cz_Moon(n + 13*i) = Cz(n);
        }
    }
    if (0<=dt && dt<=4){
        j=0;
        Mjd0 = t1;
    }
    else if(4<dt && dt<=8){
        j=1;
        Mjd0 = t1+4*j;
    }
    else if(8<dt && dt<=12){
        j=2;
        Mjd0 = t1+4*j;
    }
    else if(12<dt && dt<=16){
        j=3;
        Mjd0 = t1+4*j;
    }
    else if(16<dt && dt<=20){
        j=4;
        Mjd0 = t1+4*j;
    }
    else if(20<dt && dt<=24){
        j=5;
        Mjd0 = t1+4*j;
    }
    else if(24<dt && dt<=28){
        j=6;
        Mjd0 = t1+4*j;
    }
    else if(28<dt && dt<=32){
        j=7;
        Mjd0 = t1+4*j;
    }

    for (int n = 1; n <= 13; n++)
    {
        cx(n) = Cx_Moon(13 * j + n);
        cy(n) = Cy_Moon(13 * j + n);
        cz(n) = Cz_Moon(13 * j + n);
    }

    r_Moon = transpose(Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+4, cx,
                        cy, cz))*1e3;

    //temp = (753:11:786);
     ii = 1;
    for (int k = 753; k <= 786; k = k + 11)
    {
        temp(ii) = k;
        i++;
    }
    Matrix Cx_Sun(22);
    Matrix Cy_Sun(22);
    Matrix Cz_Sun(22);

    for (int n = 1; n <= 11; n++)
    {
        Cx_Sun(n) = PCtemp(temp(1) + n - 1);
        Cy_Sun(n) = PCtemp(temp(2) + n - 1);
        Cz_Sun(n) = PCtemp(temp(3) + n - 1);
    }

    temp = temp + 33;
    Cx = zeros(1,11);
    Cy = zeros(1,11);
    Cz = zeros(1, 11);

    for (int n = 1; n <= 11; n++)
    {
        Cx(n) = PCtemp(temp(1) + n - 1);
        Cy(n) = PCtemp(temp(2) + n - 1);
        Cz(n) = PCtemp(temp(3) + n - 1);
    }

    for (int n = 1; n <= 11; n++)
    {
        Cx_Sun(n + 11) = Cx(n);
        Cy_Sun(n + 11) = Cy(n);
        Cz_Sun(n + 11) = Cz(n);
    }

    if (0<=dt && dt<=16){
        j=0;
        Mjd0 = t1;
    }
    else if(16<dt && dt<=32){
        j=1;
        Mjd0 = t1+16*j;
    }
    for (int n = 1; n <= 11; n++)
    {
        cx(n) = Cx_Sun(11 * j + n);
        cy(n) = Cy_Sun(11 * j + n);
        cz(n) = Cz_Sun(11 * j + n);
    }
    r_Sun = transpose(Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+16, cx,
                    cy, cz))*1e3;

    //temp = (3:14:45);
     ii = 1;
    for (int k = 3; k <= 45; k = k + 14)
    {
        temp(ii) = k;
        i++;
    }

    Matrix Cx_Mercury(56);
    Matrix Cy_Mercury(56);
    Matrix Cz_Mercury(56);

    for (int n = 1; n <= 14; n++)
    {
        Cx_Mercury(n) = PCtemp(temp(1) + n - 1);
        Cy_Mercury(n) = PCtemp(temp(2) + n - 1);
        Cz_Mercury(n) = PCtemp(temp(3) + n - 1);
    }

    Cx = zeros(1, 14);
    Cy = zeros(1, 14);
    Cz = zeros(1, 14);

     for (int i = 1; i <= 3; i++)
    {
        temp = temp + 42;
        for (int n = 1; n <= 14; n++)
        {
            Cx(n) = PCtemp(temp(1) + n - 1);
            Cy(n) = PCtemp(temp(2) + n - 1);
            Cz(n) = PCtemp(temp(3) + n - 1);
        }

        for (int n = 1; n <= 14; n++)
        {
            Cx_Mercury(n + 14*i) = Cx(n);
            Cy_Mercury(n + 14*i) = Cy(n);
            Cz_Mercury(n + 14*i) = Cz(n);
        }
    }

    if (0<=dt && dt<=8){
        j=0;
        Mjd0 = t1;
    }
    else if(8<dt && dt<=16){
        j=1;
        Mjd0 = t1+8*j;
    }
    else if (16<dt && dt<=24){
        j=2;
        Mjd0 = t1+8*j;
    }
    else if(24<dt && dt<=32){
        j=3;
        Mjd0 = t1+8*j;
    }

    cx = zeros(1, 14);
    cy = zeros(1, 14);
    cz = zeros(1, 14);
    r_Mercury = transpose(Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0+8, cx,
                        cy, cz))*1e3;

    //temp = (171:10:201);
    ii = 1;
    for (int k = 171; k <= 201; k = k + 10)
    {
        temp(ii) = k;
        ii++;
    }

   Matrix Cx_Venus(20);
    Matrix Cy_Venus(20);
    Matrix Cz_Venus(20);

    for (int n = 1; n <= 10; n++)
    {
        Cx_Venus(n) = PCtemp(temp(1) + n - 1);
        Cy_Venus(n) = PCtemp(temp(2) + n - 1);
        Cz_Venus(n) = PCtemp(temp(3) + n - 1);
    }

    temp = temp + 30;

    Cx = zeros(1, 10);
    Cy = zeros(1, 10);
    Cz = zeros(1, 10);

    for (int n = 1; n <= 10; n++)
    {
        Cx(n) = PCtemp(temp(1) + n - 1);
        Cy(n) = PCtemp(temp(2) + n - 1);
        Cz(n) = PCtemp(temp(3) + n - 1);
    }

    for (int n = 1; n <= 10; n++)
    {
        Cx_Venus(n + 10) = Cx(n);
        Cy_Venus(n + 10) = Cy(n);
        Cz_Venus(n + 10) = Cz(n);
    }

    if (0<=dt && dt<=16){
        j=0;
        Mjd0 = t1;
    }
    else if(16<dt && dt<=32){
        j=1;
        Mjd0 = t1+16*j;
    }

    cx = zeros(1, 10);
    cy = zeros(1, 10);
    cz = zeros(1, 10);
    // Matrix cheb_x(10);
    // Matrix cheb_y(10);
    // Matrix cheb_z(10);
    
    for (int n = 1; n <= 10; n++)
    {
        cx(n) = Cx_Venus(10 * j + n);
        cy(n) = Cy_Venus(10 * j + n);
        cz(n) = Cz_Venus(10 * j + n);
    }

    r_Venus = transpose(Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+16, cx,
                        cy, cz))*1e3;

    // temp = (309:11:342);
    ii = 1;
    for (int k = 309; k <= 342; k = k + 11)
    {
        temp(ii) = k;
        ii++;
    }

    Matrix Cx_Mars(22);
    Matrix Cy_Mars(22);
    Matrix Cz_Mars(22);

    for (int n = 1; n <= 11; n++)
    {
        Cx_Mars(n) = PCtemp(temp(1) + n - 1);
        Cy_Mars(n) = PCtemp(temp(2) + n - 1);
        Cz_Mars(n) = PCtemp(temp(3) + n - 1);
    }

    j = 0;
    Mjd0 = t1;
    cx = zeros(1, 11);
    cy = zeros(1, 11);
    cz = zeros(1, 11);

    for (int n = 1; n <= 11; n++)
    {
        cx(n) = Cx_Mars(11 * j + n);
        cy(n) = Cy_Mars(11 * j + n);
        cz(n) = Cz_Mars(11 * j + n);
    }

    r_Mars = transpose(Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+32, cx,
                        cy, cz))*1e3;

    //temp = (342:8:366);
    ii = 1;
    for (int k = 342; k <= 366; k = k + 8)
    {
        temp(ii) = k;
        ii++;
    }

    Matrix Cx_Jupiter(16);
    Matrix Cy_Jupiter(16);
    Matrix Cz_Jupiter(16);

    for (int n = 1; n <= 8; n++)
    {
        Cx_Jupiter(n) = PCtemp(temp(1) + n - 1);
        Cy_Jupiter(n) = PCtemp(temp(2) + n - 1);
        Cz_Jupiter(n) = PCtemp(temp(3) + n - 1);
    }

    j = 0;
    Mjd0 = t1;

    cx = zeros(1, 8);
    cy = zeros(1, 8);
    cz = zeros(1, 8);
    
    for (int n = 1; n <= 8; n++)
    {
        cx(n) = Cx_Jupiter(8 * j + n);
        cy(n) = Cy_Jupiter(8 * j + n);
        cz(n) = Cz_Jupiter(8 * j + n);
    }
    r_Jupiter = transpose(Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0+32, cx,
                        cy, cz))*1e3;

    //temp = (387:6:405);
    ii = 1;
    for (int k = 366; k <= 387; k = k + 7)
    {
        temp(ii) = k;
        ii++;
    }
    Matrix Cx_Saturn(14);
    Matrix Cy_Saturn(14);
    Matrix Cz_Saturn(14);

    for (int n = 1; n <= 7; n++)
    {
        Cx_Saturn(n) = PCtemp(temp(1) + n - 1);
        Cy_Saturn(n) = PCtemp(temp(2) + n - 1);
        Cz_Saturn(n) = PCtemp(temp(3) + n - 1);
    }

    j = 0;
    Mjd0 = t1;

    cx = zeros(1, 7);
    cy = zeros(1, 7);
    cz = zeros(1, 7);
    
    for (int n = 1; n <= 7; n++)
    {
        cx(n) = Cx_Saturn(7 * j + n);
        cy(n) = Cy_Saturn(7 * j + n);
        cz(n) = Cz_Saturn(7 * j + n);
    }
    r_Saturn = transpose(Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, cx,
                        cy, cz))*1e3;

    // temp = (387 : 6 : 405);
    ii = 1;
    for (int k = 387; k <= 405; k = k + 6)
    {
        temp(ii) = k;
        ii++;
    }

    

    Matrix Cx_Uranus(12);
    Matrix Cy_Uranus(12);
    Matrix Cz_Uranus(12);

    for (int n = 1; n <= 6; n++)
    {
        Cx_Uranus(n) = PCtemp(temp(1) + n - 1);
        Cy_Uranus(n) = PCtemp(temp(2) + n - 1);
        Cz_Uranus(n) = PCtemp(temp(3) + n - 1);
    }

    j = 0;
    Mjd0 = t1;

    cx = zeros(1, 6);
    cy = zeros(1, 6);
    cz = zeros(1, 6);
    for (int n = 1; n <= 6; n++)
    {
        cx(n) = Cx_Uranus(6 * j + n);
        cy(n) = Cy_Uranus(6 * j + n);
        cz(n) = Cz_Uranus(6 * j + n);
    }

    r_Uranus = transpose( Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32, cx, cy, cz))*1e3;

    //temp = (405:6:423);
    ii = 1;
    for (int k = 405; k <= 423; k = k + 6)
    {
        temp(ii) = k;
        ii++;
    }

    Matrix Cx_Neptune(12);
    Matrix Cy_Neptune(12);
    Matrix Cz_Neptune(12);

    for (int n = 1; n <= 6; n++)
    {
        Cx_Neptune(n) = PCtemp(temp(1) + n - 1);
        Cy_Neptune(n) = PCtemp(temp(2) + n - 1);
        Cz_Neptune(n) = PCtemp(temp(3) + n - 1);
    }

    j = 0;
    Mjd0 = t1;
    
    cx = zeros(1, 6);
    cy = zeros(1, 6);
    cz = zeros(1, 6);
    
    for (int n = 1; n <= 6; n++)
    {
        cx(n) = Cx_Neptune(6 * j + n);
        cy(n) = Cy_Neptune(6 * j + n);
        cz(n) = Cz_Neptune(6 * j + n);
    }
    r_Neptune = transpose(Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, cx,
                        cy, cz))*1e3;

    //temp = (423:6:441);
     ii = 1;
    for (int k = 423; k <= 441; k = k + 6)
    {
        temp(ii) = k;
        ii++;
    }

    Matrix Cx_Pluto(12);
    Matrix Cy_Pluto(12);
    Matrix Cz_Pluto(12);

    for (int n = 1; n <= 6; n++)
    {
        Cx_Pluto(n) = PCtemp(temp(1) + n - 1);
        Cy_Pluto(n) = PCtemp(temp(2) + n - 1);
        Cz_Pluto(n) = PCtemp(temp(3) + n - 1);
    }

    j = 0;
    Mjd0 = t1;
    
    cx = zeros(1, 6);
    cy = zeros(1, 6);
    cz = zeros(1, 6);
    
    for (int n = 1; n <= 6; n++)
    {
        cx(n) = Cx_Pluto(6 * j + n);
        cy(n) = Cy_Pluto(6 * j + n);
        cz(n) = Cz_Pluto(6 * j + n);
    }
    r_Pluto = transpose(Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, cx,
                        cy, cz))*1e3;

    /*temp = (819:10:839);
    Cx_Nutations = PCtemp(temp(1):temp(2)-1);
    Cy_Nutations = PCtemp(temp(2):temp(3)-1);
    for i=1:3
        temp = temp+20;
        Cx = PCtemp(temp(1):temp(2)-1);
        Cy = PCtemp(temp(2):temp(3)-1);
        Cx_Nutations = [Cx_Nutations,Cx];
        Cy_Nutations = [Cy_Nutations,Cy];
    }
    if (0<=dt && dt<=8)
        j=0;
        Mjd0 = t1;
    elseif(8<dt && dt<=16)
        j=1;
        Mjd0 = t1+8*j;
    elseif (16<dt && dt<=24)
        j=2;
        Mjd0 = t1+8*j;
    elseif(24<dt && dt<=32)
        j=3;
        Mjd0 = t1+8*j;
    }
    Nutations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, Cx_Nutations(10*j+1:10*j+10),...
                    Cy_Nutations(10*j+1:10*j+10),zeros(10,1))';

    temp = (899:10:929);
    Cx_Librations = PCtemp(temp(1):temp(2)-1);
    Cy_Librations = PCtemp(temp(2):temp(3)-1);
    Cz_Librations = PCtemp(temp(3):temp(4)-1);
    for i=1:3
        temp = temp+30;
        Cx = PCtemp(temp(1):temp(2)-1);
        Cy = PCtemp(temp(2):temp(3)-1);
        Cz = PCtemp(temp(3):temp(4)-1);
        Cx_Librations = [Cx_Librations,Cx];
        Cy_Librations = [Cy_Librations,Cy];
        Cz_Librations = [Cz_Librations,Cz];    
    }
    if (0<=dt && dt<=8)
        j=0;
        Mjd0 = t1;
    elseif(8<dt && dt<=16)
        j=1;
        Mjd0 = t1+8*j;
    elseif (16<dt && dt<=24)
        j=2;
        Mjd0 = t1+8*j;
    elseif(24<dt && dt<=32)
        j=3;
        Mjd0 = t1+8*j;
    }
    Librations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, Cx_Librations(10*j+1:10*j+10),...
                        Cy_Librations(10*j+1:10*j+10), Cz_Librations(10*j+1:10*j+10))';*/
    double EMRAT = 81.30056907419062; // DE430
    double EMRAT1 = 1/(1+EMRAT);
    r_Earth = r_Earth-r_Moon*EMRAT1;
    r_Mercury = r_Mercury-r_Earth;
    r_Venus =  r_Venus-r_Earth;
    r_Mars =  r_Mars-r_Earth;
    r_Jupiter =  r_Jupiter-r_Earth;
    r_Saturn =  r_Saturn-r_Earth;
    r_Uranus =  r_Uranus-r_Earth;
    r_Neptune =  r_Neptune-r_Earth;
    r_Pluto =  r_Pluto-r_Earth;
    r_Sun =  r_Sun-r_Earth;
}
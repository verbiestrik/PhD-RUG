function [r,theta,h,vr,alpha1,alpha2,alpha3,vtheta,gamma1,gamma2,gamma3] = readDataASWME(name,type)
if(strcmp(type,'ref'))
    matrix = load(strcat(name,'.csv'));
    r = matrix(:,1);
    theta = matrix(:,2);
    h = matrix(:,3);
    hvr =  matrix(:,4);
    halpha1 =  matrix(:,5);
    halpha2 =  matrix(:,6);
    halpha3 =  matrix(:,7);
    hvtheta =  matrix(:,8);
    hgamma1 =  matrix(:,9);
    hgamma2 =  matrix(:,10);
    hgamma3 =  matrix(:,11);

    vr = hvr;
    alpha1 =  halpha1;
    alpha2 =  halpha2;
    alpha3 =  halpha3;
    vtheta =  hvtheta;
    gamma1 =  hgamma1;
    gamma2 =  hgamma2;
    gamma3 =  hgamma3;
    
%     vr = hvr ./ h;
%     alpha1 =  halpha1 ./ h;
%     alpha2 =  halpha2 ./ h;
%     alpha3 =  halpha3 ./ h;
%     vtheta =  hvtheta ./ h;
%     gamma1 =  hgamma1 ./ h;
%     gamma2 =  hgamma2 ./ h;
%     gamma3 =  hgamma3 ./ h;

    alpha1 =  3. * alpha1;
    alpha2 =  5. * alpha2;
    alpha3 =  7. * alpha3;

    gamma1 =  3. * gamma1;
    gamma2 =  5. * gamma2;
    gamma3 =  7. * gamma3;

end
if(strcmp(type,'swe'))
    matrix = load(strcat(name,'.csv'));
    r = matrix(:,1);
    theta = matrix(:,2);
    h = matrix(:,3);
    hvr =  matrix(:,4);
    hvtheta = matrix(:,5);
    alpha1 =  0*matrix(:,5);
    alpha2 =  0*matrix(:,5);
    alpha3 =  0*matrix(:,5);
    gamma1 =  0*matrix(:,5);
    gamma2 =  0*matrix(:,5);
    gamma3 =  0*matrix(:,5);

    vtheta = hvtheta ./ h;
    vr = hvr ./ h;
end
if(strcmp(type,'swme'))
    matrix = load(strcat(name,'.csv'));
    r = matrix(:,1);
    theta = matrix(:,2);
    h = matrix(:,3);
    hvr =  matrix(:,4);
    if(size(matrix,2) == 7)
        halpha1 =  matrix(:,5);
        hvtheta =  matrix(:,6);
        hgamma1 =  matrix(:,7);

        halpha2 =  0 * matrix(:,7);
        halpha3 =  0 * matrix(:,7);
        hgamma2 =  0 * matrix(:,7);
        hgamma3 =  0 * matrix(:,7);
    end
    if(size(matrix,2) == 9)
        halpha1 =  matrix(:,5);
        halpha2 =  matrix(:,6);
        hvtheta =  matrix(:,7);
        hgamma1 =  matrix(:,8);
        hgamma2 =  matrix(:,9);

        halpha3 =  0 * matrix(:,9);
        hgamma3 =  0 * matrix(:,9);
    end
    if(size(matrix,2) == 11)
        halpha1 =  matrix(:,5);
        halpha2 =  matrix(:,6);
        halpha3 =  matrix(:,7);
        hvtheta =  matrix(:,8);
        hgamma1 =  matrix(:,9);
        hgamma2 =  matrix(:,10);
        hgamma3 =  matrix(:,11);
    end
    if(size(matrix,2) == 13)
        halpha1 =  matrix(:,5);
        halpha2 =  matrix(:,6);
        halpha3 =  matrix(:,7);
        hvtheta =  matrix(:,9);
        hgamma1 =  matrix(:,10);
        hgamma2 =  matrix(:,11);
        hgamma3 =  matrix(:,12);
    end
    if(size(matrix,2) == 15)
        halpha1 =  matrix(:,5);
        halpha2 =  matrix(:,6);
        halpha3 =  matrix(:,7);
        hvtheta =  matrix(:,10);
        hgamma1 =  matrix(:,11);
        hgamma2 =  matrix(:,12);
        hgamma3 =  matrix(:,13);
    end
    vr = hvr ./ h;
    alpha1 =  halpha1 ./ h;
    alpha2 =  halpha2 ./ h;
    alpha3 =  halpha3 ./ h;
    vtheta =  hvtheta ./ h;
    gamma1 =  hgamma1 ./ h;
    gamma2 =  hgamma2 ./ h;
    gamma3 =  hgamma3 ./ h;
%    test = u + alpha1 + alpha2; % for EQ 3 plots
%     test = alpha2; % for EQ 3 plots
%    alpha2 = test;
    

%     alpha1 =  3 * alpha1;
%     alpha2 =  5 * alpha2;
end
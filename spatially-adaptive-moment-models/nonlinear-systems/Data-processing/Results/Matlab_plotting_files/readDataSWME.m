function [x,h,u,alpha1,alpha2] = readDataSWME(name,type)
if(strcmp(type,'ref'))
    matrix = load(strcat(name,'.csv'));
    x = matrix(:,1);
    h = matrix(:,2);
    hu =  matrix(:,3);
    halpha1 =  matrix(:,4);
    halpha2 =  matrix(:,5);

    u = hu;
    alpha1 =  halpha1;
    alpha2 =  halpha2;
    
%     u = hu ./ h;
%     alpha1 =  halpha1 ./ h;
%     alpha2 =  halpha2 ./ h;

    alpha1 =  3 * alpha1;
    alpha2 =  5 * alpha2;
end
if(strcmp(type,'swe'))
    matrix = load(strcat(name,'.csv'));
    x = matrix(:,1);
    h = matrix(:,2);
    hu =  matrix(:,3);
    alpha1 =  0*matrix(:,3);
    alpha2 =  0*matrix(:,3);

    u = hu ./ h;
end
if(strcmp(type,'swme'))
    matrix = load(strcat(name,'.csv'));
    x = matrix(:,1);
    h = matrix(:,2);
    hu =  matrix(:,3);
    halpha1 =  matrix(:,4);
    if(size(matrix,2)>4)
        halpha2 =  matrix(:,5);
    else
        halpha2 = 0 * matrix(:,4);
    end
    u = hu ./ h;
    alpha1 =  halpha1 ./ h;
    alpha2 =  halpha2 ./ h;
    test = u + alpha1 + alpha2; % for EQ 3 plots
%     test = alpha2; % for EQ 3 plots
    alpha2 = test;
    

%     alpha1 =  3 * alpha1;
%     alpha2 =  5 * alpha2;
end
if(strcmp(type,'swmep'))
    matrix = load(strcat(name,'.csv'));
    x = matrix(:,1);
    h = matrix(:,2);
    u =  matrix(:,3);
    alpha1 =  matrix(:,4);
    if(size(matrix,2)>4)
        alpha2 =  matrix(:,5);
    else
        alpha2 = 0 * matrix(:,4);
    end
end
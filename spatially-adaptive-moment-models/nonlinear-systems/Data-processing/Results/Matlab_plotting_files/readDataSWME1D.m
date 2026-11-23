function [x,h,u,alpha1,alpha2,alpha3,alpha4,alpha5] = readDataSWME(name,type)
if(strcmp(type,'ref'))
    matrix = load(strcat(name,'.csv'));
    x = matrix(:,1);
    h = matrix(:,2);
    hu =  matrix(:,3);
    halpha1 =  matrix(:,4);
    halpha2 =  matrix(:,5);
    halpha3 =  matrix(:,6);
    halpha4 =  matrix(:,7);
    halpha5 =  matrix(:,8);

    % If the outputs of the numerical reference solver are the primitive variables:
    u = hu;
    alpha1 =  halpha1;
    alpha2 =  halpha2;
    alpha3 =  halpha3;
    alpha4 =  halpha4;
    alpha5 =  halpha5;
    
    % If the outputs of the numerical reference solver are not the primitive variables:
%     u = hu ./ h;
%     alpha1 =  halpha1 ./ h;
%     alpha2 =  halpha2 ./ h;
%     alpha3 =  halpha3 ./ h;
%     alpha4 =  halpha4 ./ h;
%     alpha5 =  halpha5 ./ h;

    alpha1 =  3 * alpha1;
    alpha2 =  5 * alpha2;
    alpha3 =  7 * alpha3;
    alpha4 =  9 * alpha4;
    alpha5 =  11 * alpha5;
end
if(strcmp(type,'swe'))
    matrix = load(strcat(name,'.csv'));
    x = matrix(:,1);
    h = matrix(:,2);
    hu =  matrix(:,3);
    alpha1 =  0*matrix(:,3);
    alpha2 =  0*matrix(:,3);
    alpha3 =  0*matrix(:,3);
    alpha4 =  0*matrix(:,3);
    alpha5 =  0*matrix(:,3);
    
    % If the outputs of the numerical solver are the primitive variables:
    u = hu;

    % If the outputs of the numerical solver are not the primitive variables:
    %u = hu ./ h;
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
    if(size(matrix,2)>5)
        halpha3 =  matrix(:,6);
    else
        halpha3 = 0 * matrix(:,4);
    end
    if(size(matrix,2)>6)
        halpha4 =  matrix(:,7);
    else
        halpha4 = 0 * matrix(:,4);
    end
    if(size(matrix,2)>7)
        halpha5 =  matrix(:,8);
    else
        halpha5 = 0 * matrix(:,4);
    end

    % If the outputs of the numerical solver are the primitive variables:
    u = hu;
    alpha1 =  halpha1;
    alpha2 =  halpha2;
    alpha3 =  halpha3;
    alpha4 =  halpha4;
    alpha5 =  halpha5;
    
    % If the outputs of the numerical solver are not the primitive variables:
    %u = hu./ h;
    %alpha1 =  halpha1./ h;
    %alpha2 =  halpha2./ h;
    %alpha3 =  halpha3./ h;
    %alpha4 =  halpha4./ h;
    %alpha5 =  halpha5./ h;

end
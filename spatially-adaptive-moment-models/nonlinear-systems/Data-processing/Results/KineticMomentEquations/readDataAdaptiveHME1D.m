function [x,rho,u,theta,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,moments] = readDataAdaptiveHME1D(name,type)
if(strcmp(type,'hme'))
    matrix = load(strcat(name,'.csv'));
    x = matrix(:,1);
    rho = matrix(:,2);
    u =  matrix(:,3);
    theta =  matrix(:,4);
    if(size(matrix,2)>4)
        f3 =  matrix(:,5);
    else
        f3 = 0 * matrix(:,4);
    end
    if(size(matrix,2)>5)
        f4 =  matrix(:,6);
    else
        f4 = 0 * matrix(:,4);
    end
    if(size(matrix,2)>6)
        f5 =  matrix(:,7);
    else
        f5 = 0 * matrix(:,4);
    end
    if(size(matrix,2)>7)
        f6 =  matrix(:,8);
    else
        f6 = 0 * matrix(:,4);
    end
    if(size(matrix,2)>8)
        f7 =  matrix(:,9);
    else
        f7 = 0 * matrix(:,4);
    end
    if(size(matrix,2)>9)
        f8 =  matrix(:,10);
    else
        f8 = 0 * matrix(:,4);
    end
    if(size(matrix,2)>10)
        f9 =  matrix(:,11);
    else
        f9 = 0 * matrix(:,4);
    end
    if(size(matrix,2)>11)
        f10 =  matrix(:,12);
    else
        f10 = 0 * matrix(:,4);
    end
    if(size(matrix,2)>12)
        f11 =  matrix(:,13);
    else
        f11 = 0 * matrix(:,4);
    end
    if(size(matrix,2)>13)
        f12 =  matrix(:,14);
    else
        f12 = 0 * matrix(:,4);
    end
    
    if size(matrix,2) == 4
        moments = 2*ones(size(x));
    elseif size(matrix,2) == 5
        moments = 3*ones(size(x));
    elseif size(matrix,2) == 6
        moments = 4*ones(size(x));
    elseif size(matrix,2) == 7
        moments = 5*ones(size(x));
    elseif size(matrix,2) == 8
        moments = 6*ones(size(x));
    elseif size(matrix,2) == 9
        moments = 7*ones(size(x));
    elseif size(matrix,2) == 10
        moments = 8*ones(size(x));
    elseif size(matrix,2) == 11
        moments = 9*ones(size(x));
    elseif size(matrix,2) == 12
        moments = 10*ones(size(x));
    elseif size(matrix,2) == 13
        moments = 11*ones(size(x));
    elseif size(matrix,2) == 14
        moments = 12*ones(size(x));
    end
end

if(strcmp(type,'hmeAdaptive'))
    matrix = load(strcat(name,'.csv'));
    x = matrix(:,1);
    rho = matrix(:,2);
    u =  matrix(:,3);
    theta =  matrix(:,4);
    if(size(matrix,2)>4)
        f3 =  matrix(:,5);
    else
        f3 = 0 * matrix(:,4);
    end
    if(size(matrix,2)>5)
        f4 =  matrix(:,6);
    else
        f4 = 0 * matrix(:,4);
    end
    if(size(matrix,2)>6)
        f5 =  matrix(:,7);
    else
        f5 = 0 * matrix(:,4);
    end
    if(size(matrix,2)>7)
        f6 =  matrix(:,8);
    else
        f6 = 0 * matrix(:,4);
    end
    if(size(matrix,2)>8)
        f7 =  matrix(:,9);
    else
        f7 = 0 * matrix(:,4);
    end
    if(size(matrix,2)>9)
        f8 =  matrix(:,10);
    else
        f8 = 0 * matrix(:,4);
    end
    if(size(matrix,2)>10)
        f9 =  matrix(:,11);
    else
        f9 = 0 * matrix(:,4);
    end
    if(size(matrix,2)>11)
        f10 =  matrix(:,12);
    else
        f10 = 0 * matrix(:,4);
    end
    if(size(matrix,2)>12)
        f11 =  matrix(:,13);
    else
        f11 = 0 * matrix(:,4);
    end
    if(size(matrix,2)>13)
        f12 =  matrix(:,14);
    else
        f12 = 0 * matrix(:,4);
    end
    
    moments = matrix(:,end);
end
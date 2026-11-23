function [data_matrix] = readDataHME_matrixForm(name)
matrix = load(strcat(name,'.csv'));
x = matrix(:,1);
data_matrix = zeros(length(x),12);
data_matrix(:,1) = matrix(:,1);
data_matrix(:,2) = matrix(:,2);
data_matrix(:,3) =  matrix(:,3);
data_matrix(:,4) =  matrix(:,4);
if(size(matrix,2)>4)
    data_matrix(:,5) =  matrix(:,5);
else
    data_matrix(:,5) = 0 * matrix(:,4);
end
if(size(matrix,2)>5)
    data_matrix(:,6) =  matrix(:,6);
else
    data_matrix(:,6) = 0 * matrix(:,4);
end
if(size(matrix,2)>6)
    data_matrix(:,7) =  matrix(:,7);
else
    data_matrix(:,7) = 0 * matrix(:,4);
end
if(size(matrix,2)>7)
    data_matrix(:,8) =  matrix(:,8);
else
    data_matrix(:,8) = 0 * matrix(:,4);
end
if(size(matrix,2)>8)
    data_matrix(:,9) =  matrix(:,9);
else
    data_matrix(:,9) = 0 * matrix(:,4);
end
if(size(matrix,2)>9)
    data_matrix(:,10) =  matrix(:,10);
else
    data_matrix(:,10) = 0 * matrix(:,4);
end
if(size(matrix,2)>10)
    data_matrix(:,11) =  matrix(:,11);
else
    data_matrix(:,11) = 0 * matrix(:,4);
end
if(size(matrix,2)>11)
    data_matrix(:,12) =  matrix(:,12);
else
    data_matrix(:,12) = 0 * matrix(:,4);
end
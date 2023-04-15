function [nb_num, nb_face] = findnb(eNum, mypath)
% findnb to find a list of neighbors of current element and the
% corresponding faces
% Input: eNum - current element number
% Input: mypath - build path
% Output: nb_num - list of nb elements number
% Output: nb_face - nb faces


dis = mypath(1:eNum-1,:) - mypath(eNum,:);
dis_abs = abs(dis);
nl = (sum(dis_abs,2) == 1); % logical vector to select neighbours
tmpn = 1:eNum-1;
nb_num = tmpn(nl); % List of neighbours
% since only considering immediate neighbours here, so the sum of abs(dis)
% must equal to 1 to qualify as a neighbour.
nl_dis = dis(nl,:);

%dis = sqrt(sum(dis.^2,2))
mat1 = [1 0 0; 0 1 0; 0 0 1];
mat1 = [mat1; -mat1];

c = logical(mat1*nl_dis' == 1);
tmp_contact = repmat((1:6)',1,size(c,2));
nb_face = tmp_contact(c)'; % List of contact face number

end



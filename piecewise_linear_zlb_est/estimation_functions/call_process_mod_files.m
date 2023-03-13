% modnam_00 = 'BNK_model'; % base model (constraint 1 and 2 below don't bind)
% modnam_10 = 'BNK_model_zlb';  % first constraint is true
% % modnam_10 = 'BNK_model';  % first constraint is true
% modnam_01 = 'BNK_model';  % second constraint is true
% modnam_11 = 'BNK_model'; % both constraints bind
% constraint1       = 'r_t  < r_zlb'; 
% constraint_relax1 = 'rnot > r_zlb';
% % constraint2       = 'r_t  < r_zlb'; 
% % constraint_relax2 = 'rnot > r_zlb';
% constraint2       = 'pi_t  < 10'; 
% constraint_relax2 = 'pi_t > 10';

warning('off','all')
warning

curb_retrench =0;
maxiter = 10;

eval(['dynare ',modnam_00,' nolog noclearall'])

wishlist_ = M_.endo_names;
nwishes_ = M_.endo_nbr;

obs_list_withrnot = union(obs_list,'data_rnot','rows');
obs_list_nor = setdiff(obs_list_withrnot,{'data_rnot';'data_r'},'rows');
obs_list_r = intersect(obs_list_withrnot,{'data_rnot';'data_r'},'rows');

[~, i1, ~]=intersect(wishlist_,obs_list_nor,'rows');
[~, i2, ~]=intersect(wishlist_,obs_list_r,'rows');
disp('create zdata script')
fid = fopen('eval_zdata_script.m','wt');
for i_indx_=i1'
  fprintf(fid,[deblank(wishlist_(i_indx_,:)),'_l=zdatal(:,' num2str(i_indx_) ');\n']);
  fprintf(fid,[deblank(wishlist_(i_indx_,:)),'_p=zdatap(:,' num2str(i_indx_) ');\n']);
  fprintf(fid,[deblank(wishlist_(i_indx_,:)),'_ss=zdatass(' num2str(i_indx_) ');\n']);
end
for i_indx_=i2'
  fprintf(fid,[deblank(wishlist_(i_indx_,:)),'_l=exp(zdatal(:,' num2str(i_indx_) '))-1;\n']);
  fprintf(fid,[deblank(wishlist_(i_indx_,:)),'_p=exp(zdatap(:,' num2str(i_indx_) '))-1;\n']);
  fprintf(fid,[deblank(wishlist_(i_indx_,:)),'_ss=zdatass(' num2str(i_indx_) ');\n']);
end
fclose(fid);


fid = fopen('eval_endo_names.m','wt');
for i_indx_=1:M_.endo_nbr
  fprintf(fid,[deblank(M_.endo_names(i_indx_,:)),'_ss=oo00_.dr.ys(' num2str(i_indx_) ');\n']);
end
fclose(fid);

fid = fopen('eval_param.m','wt');
for i_indx_=1:size(M_.param_names,1);
  fprintf(fid,[deblank(M_.param_names(i_indx_,:)),'=M00_.params(' num2str(i_indx_) ');\n']);
end
fclose(fid);




%---------------------------------------
% Create all the matrices used by OccBin
% processes the constraint so as to uppend a suffix to each
%---------------------------------------
solve_two_constraints_firstcall(modnam_00,modnam_10,modnam_01,modnam_11);

[constraint1_difference iendo1]= process_constraint_with_tokens(constraint1,'_difference',M00_.endo_names,0);
[constraint_relax1_difference iendo2]= process_constraint_with_tokens(constraint_relax1,'_difference',M00_.endo_names,0);
[constraint2_difference iendo3]= process_constraint_with_tokens(constraint2,'_difference',M00_.endo_names,0);
[constraint_relax2_difference iendo4]= process_constraint_with_tokens(constraint_relax2,'_difference',M00_.endo_names,0);

iendo_constraint = union(union(union(iendo1,iendo2),iendo3),iendo4);

fid = fopen('eval_difference_script.m','wt');
for i_indx_=iendo_constraint
	fprintf(fid,[deblank(wishlist_(i_indx_,:)),'_difference=zdatalinear_(:,' num2str(i_indx_) ');\n']);
end
fclose(fid);


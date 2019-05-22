%% Section 1: Run this section first before proceeding
clear;
clc;
close all;
E0=200*1e3;%accelerating voltage in V
P_rotated=imread('centre from bottom 11.tif');
[totrow,totcol]=size(P_rotated);
Q_rotated=zeros(totrow+400,totcol+400);% padding to avoid errors
Q_rotated(1:totrow,1:totcol)=P_rotated(1:totrow,1:totcol);
P_rotated=Q_rotated;
load('I400_5mrad.mat');
% Rough knowledge about the In composition, X_In, you can choose from X_In=0.00, 0.25, 0.50, 0.75, 1.00
%setting the scalebar
X_In=0.50;% roughly estimated composition
comp_number=3;% put 1 when X_In=0   put 2 if X_In=0.25   put 3 if X_In=0.50  put 4 if X_In=0.75 put 5 if X_In=1.00 
scalebar_initial=36.9;%pixels/(1/nm) keep 3 digits at most @ 300mm F20 this is 70.9 @ Titan this is 37.1/37.2 @Titan Aug22 35.8 @Oct 9(300kV) 33 Oct 18-19 (200kV)37.4/38.9 newly measured 

%% Section 2.1: find the maximum point from the pattern to help locating the central disc
[maxValue, linearIndexesOfMaxes] = max(P_rotated(:));
[rowsOfMaxes, colsOfMaxes] = find(P_rotated == maxValue);

scalebar=scalebar_initial*10; %pixels/(1/A)
g400length=g400_length(X_In); %length of g400 in 1/A
g400length=g400length*scalebar; %length of g400 in # of pixels
g200length=0.5*g400length; 
%initially, we draw a horizontal line across the point of maximum
initial_row=P_rotated(rowsOfMaxes,round(colsOfMaxes-2*g200length):round(colsOfMaxes+2*g200length));
figure
plot(initial_row);
initial_col=P_rotated(round(rowsOfMaxes-2*g200length):round(rowsOfMaxes+2*g200length),colsOfMaxes);
figure
plot(initial_col);
[~,size_initial_row]=size(initial_row);
d_initial_row=zeros(1,size_initial_row);

for n=1:size_initial_row-1
     d_initial_row(n+1)=initial_row(n+1)-initial_row(n);
end
figure
plot(d_initial_row);

%find the maximum and return the position of the peak
[max_num, ~]=max(d_initial_row(:));
[X_leftedge]=ind2sub(size(d_initial_row),find(d_initial_row==max_num));
%find the minimum and return the position of the valley
[min_num, ~]=min(d_initial_row(:));
[X_rightedge]=ind2sub(size(d_initial_row),find(d_initial_row==min_num));
X_centre=0.5*(X_rightedge+X_leftedge);

%% Section 2.2 Find the centre of the 000 disc and find the diametre of the disc
[size_initial_col,~]=size(initial_col);
d_initial_col=zeros(size_initial_col,1);

for n=1:size_initial_col-1
     d_initial_col(n+1)=initial_col(n+1)-initial_col(n);
end

figure
plot(d_initial_col);
title('column edge')
%find the maximum and return the position of the peak
[max_num2, ~]=max(d_initial_col(:));
[Y_topedge]=ind2sub(size(d_initial_col),find(d_initial_col==max_num2));
%find the minimum and return the position of the valley
[min_num2, ~]=min(d_initial_col(:));
[Y_bottomedge]=ind2sub(size(d_initial_col),find(d_initial_col==min_num2));
Y_centre=0.5*(Y_bottomedge+Y_topedge);
%? Not sure if the coordinate has an error of 1 pixel
col_centre=round(X_centre+colsOfMaxes-2*g200length-1);%avoid factional pixel; the cordinate of the centre of the 000 disc
row_centre=round(Y_centre+rowsOfMaxes-2*g200length-1);
% we draw a vertical line and horizontal line across the centre of the 000
% disc

centre_row=P_rotated(row_centre,col_centre-round(2*g200length):col_centre+round(2*g200length));
figure
plot(centre_row);
title('centre row');
[~,size_centre_row]=size(centre_row);
d_centre_row=zeros(1,size_centre_row);

for n=1:size_centre_row-1
     d_centre_row(n+1)=centre_row(n+1)-centre_row(n);
end
figure
plot(d_centre_row);
title('d_centre_row');
%find the maximum and return the position of the peak
[max_num3, ~]=max(d_centre_row(:));
[X_centrerow_leftedge]=ind2sub(size(d_centre_row),find(d_centre_row==max_num3));
%find the minimum and return the position of the valley
[min_num3, ~]=min(d_centre_row(:));
[X_centrerow_rightedge]=ind2sub(size(d_centre_row),find(d_centre_row==min_num3));
diametre=X_centrerow_rightedge-X_centrerow_leftedge;
diametre_initial=diametre;
%% Section 3.1 capture the -400 and 400 disc
X_tl_left400=round(col_centre-g400length-0.50*diametre);
Y_tl_left400=round(row_centre-0.50*diametre);
X_br_left400=round(col_centre-g400length+0.50*diametre);
Y_br_left400=round(row_centre+0.50*diametre);
left400=P_rotated(Y_tl_left400:Y_br_left400,X_tl_left400:X_br_left400);
left400_upper_bg=P_rotated(round(Y_tl_left400-0.50*diametre):Y_tl_left400,X_tl_left400:X_br_left400);
left400_lower_bg=P_rotated(Y_br_left400:round(Y_br_left400+0.5*diametre),X_tl_left400:X_br_left400);

Displayed_left400=log10(1+0.005*abs(left400));
figure
imshow(Displayed_left400);
title('-400')
figure
left400_profile=sum(left400,1)-sum(left400_upper_bg,1)-sum(left400_lower_bg,1);

plot(left400_profile);
title('profile of -400')

X_tl_right400=round(col_centre+g400length-0.50*diametre);
Y_tl_right400=round(row_centre-0.50*diametre);
X_br_right400=round(col_centre+g400length+0.50*diametre);
Y_br_right400=round(row_centre+0.50*diametre);

right400=P_rotated(Y_tl_right400:Y_br_right400,X_tl_right400:X_br_right400);
right400_upper_bg=P_rotated(round(Y_tl_right400-0.50*diametre):Y_tl_right400,X_tl_right400:X_br_right400);
right400_lower_bg=P_rotated(Y_br_right400:round(Y_br_right400+0.5*diametre),X_tl_right400:X_br_right400);
Displayed_right400=log10(1+0.005*abs(right400));

figure
imshow(Displayed_right400);
title('400');

right400_profile=sum(right400,1)-sum(right400_upper_bg,1)-sum(right400_lower_bg,1);
figure
plot(right400_profile);
title('profile of 400');
figure
%% Section 4.1 now draw a rectangle that covers the 400 kicuchi bands and record the top-left and bottom-right corners of the
%rectangle
x_br=round(col_centre+g200length+0.5*diametre);%br means bottom left
y_br=round(row_centre-g400length);
x_tl=round(x_br-g400length-diametre);%tl means top left
y_tl=round(y_br-g400length);

row_max=round(row_centre-g400length-diametre);
row_min=round(row_max-diametre-2*g400length);
  if row_min<0
    row_min=1;
  end

x_max=x_br-x_tl+1;
y_max=y_br-y_tl+1;
crop_kikuchi=zeros(y_max,x_max);

for x=1:x_max
    for y=1:y_max
        crop_kikuchi(y,x)=P_rotated(y+y_tl-1,x+x_tl-1);
    end
end
Displayed_crop_kikchi=log10(1+0.005*abs(crop_kikuchi));
figure
imshow(Displayed_crop_kikchi);

%% Section 4.2: Take the profile of the summed kikuchi band
Kikuchi_profile=sum(crop_kikuchi);
Kikuchi_profile=smooth(Kikuchi_profile);
figure
plot(1:x_max,Kikuchi_profile);
smoothed_Kikuchi=smooth(Kikuchi_profile);
%smoothed_Kikuchi=smoothed_Kikuchi';
[size_smoothed,~]=size(smoothed_Kikuchi);
figure
plot(1:size_smoothed,smoothed_Kikuchi);
%% Section 4.3: Take the first-order derivative
first_order=zeros(size_smoothed,1);
 first_order(1)=0;
for n=1:size_smoothed-1
     first_order(n+1)=smoothed_Kikuchi(n+1)-smoothed_Kikuchi(n);
end
figure
plot(1:size_smoothed,first_order);
first_order=smooth(smooth(smooth(first_order)));
figure
plot(1:size_smoothed,first_order);
%find the maximum and return the position of the peak
[max_num, ~]=max(first_order(:));
[X_L_Kikuchi]=ind2sub(size(first_order),find(first_order==max_num));
%find the minimum and return the position of the valley
[min_num, ~]=min(first_order(:));
[X_R_Kikuchi]=ind2sub(size(first_order),find(first_order==min_num));
%% Section 4.4: Find the edges of the 400 Kikuchi band

%start from here
threshold=0.5*min(abs(max_num),abs(min_num));
Binary_profile = imbinarize(abs(first_order),threshold);
figure
plot(Binary_profile);

profile_bin=2*Binary_profile;

locs_edge = find(profile_bin);

imax=numel(locs_edge);
edge_candidate_matrix=zeros(imax,4);
col_candidate=1;

edgesize=zeros(1,4);
% edge_size_i=zeros(1,4);% an array of the widths of the candidate edges/flattened peaks

    for i=1:imax-1
       if locs_edge(i+1)-locs_edge(i)<=2
         edge_candidate_matrix(i,col_candidate)=locs_edge(i);
       else
           col_candidate=col_candidate+1;
       end
    end
   number_of_pixel_peredge=zeros(1,4);
    for edge_number=1:4
    number_of_pixel_peredge(edge_number)=nnz(edge_candidate_matrix(:,edge_number));
    end
    [temp_rank,original_pos]=sort(number_of_pixel_peredge,'descend');
two_edges=original_pos(1:2);%take the first two candiate edge_numbers
if two_edges(1)>two_edges(2) % put the left edge on the left
    temp=two_edges(1);
    two_edges(1)=two_edges(2);
    two_edges(2)=temp;
end
left_edge_col=two_edges(1);
right_edge_col=two_edges(2);
X_L_Kikuchi=sum(edge_candidate_matrix(:,left_edge_col))/nnz(edge_candidate_matrix(:,left_edge_col));
X_R_Kikuchi=sum(edge_candidate_matrix(:,right_edge_col))/nnz(edge_candidate_matrix(:,right_edge_col));

Kikuchi_width=X_R_Kikuchi-X_L_Kikuchi; % in # of pixels
centre_Kikuchi=0.5*(X_R_Kikuchi+X_L_Kikuchi);

centre_Kikuchi=centre_Kikuchi+x_tl;
% the deviation from mirror symmetry
delta_X=col_centre-centre_Kikuchi; % in # of pixels
lambda=12.2643/sqrt(E0*(1+0.978476*1e-6*E0)); % the relativiscally corrected wavlength

scalebar=scalebar_initial*10; %pixels/(1/A)
scalebar_mrad=scalebar*1e-3/lambda; %pixels/mrad
scalebar_mrad=round(scalebar_mrad); %pixels(integer)/mrad
delta_X_mrad=delta_X/scalebar_mrad; % in mrad
radius_mrad=0.5*diametre/scalebar_mrad; % in mrad
%% Section 5: extract the precalculated thickness profiles and save them within the diametre of the disc

scalebar_cal=17; % Assuming that I400_5mrad.mat is loaded, so the scalebar is 17 pixels/mrad

cross_correlation=zeros(200,41);
difference=1000*ones(200,41);
likelihood=zeros(200,41);
% difference=zeros(200,41);
diff_pks=1000*ones(200,41);
for i=1:200
    
[~,sizeleft400]=size(left400_profile);%experiment
[~,sizeright400]=size(right400_profile);


 for dx=1:41

  I400_left_new=I400_left(:,round(dx-11+86+delta_X_mrad*scalebar_cal-radius_mrad*scalebar_cal):round(dx-21+86+delta_X_mrad*scalebar_cal+radius_mrad*scalebar_cal),comp_number); % calculation
  I400_right_new=I400_right(:,round(dx-11+86+delta_X_mrad*scalebar_cal-radius_mrad*scalebar_cal):round(dx-21+86+delta_X_mrad*scalebar_cal+radius_mrad*scalebar_cal),comp_number);

%% Section 6: resize the experimental I-400 and I400 profile


[~,sizeleft400_cal,~]=size(I400_left_new); %calculation
[~,sizeright400_cal,~]=size(I400_right_new);

left400_profile_new=imresize(left400_profile,sizeleft400_cal/sizeleft400); % resize the experimental profile to the same scale as the calculated one
right400_profile_new=imresize(right400_profile,sizeright400_cal/sizeright400);
left400_exp=left400_profile_new(1,:)/(maxValue*diametre_initial); % divided by the maximum intensity for normalization and take the first row (don't know why there were two rows in lef_profile_new)
right400_exp=right400_profile_new(1,:)/(maxValue*diametre_initial); 

%% cross-correlation between experiments and calculation

left400_exp=left400_exp(1,35:sizeleft400-20); % avoid the edge
right400_exp=right400_exp(1,15:sizeright400-40);

left400_cal=I400_left_new(i,35:sizeleft400-20)/171;
right400_cal=I400_right_new(i,15:sizeright400-40)/171;

    
% cross_correlation_tensor(i,dx)=sumsqr(xcorr(left400_cal,left400_exp))+sumsqr(xcorr(right400_cal,right400_exp));
% difference(i,dx)=sumsqr(left400_cal-left400_exp)+sumsqr(right400_cal-right400_exp);

[ltpks,ltpks_locs]=findpeaks(left400_exp);
[ltvls,ltvls_locs]=findpeaks(-left400_exp);
[rtpks,rtpks_locs]=findpeaks(right400_exp);
[rtvls,rtvls_locs]=findpeaks(-right400_exp);
[~,size_ltpks]=size(ltpks_locs);
[~,size_ltvls]=size(ltvls_locs);
[~,size_rtpks]=size(rtpks_locs);
[~,size_rtvls]=size(rtvls_locs);


[ltpks_cal,ltpks_locs_cal]=findpeaks(left400_cal);
[ltvls_cal,ltvls_locs_cal]=findpeaks(-left400_cal);
[rtpks_cal,rtpks_locs_cal]=findpeaks(right400_cal);
[rtvls_cal,rtvls_locs_cal]=findpeaks(-right400_cal);

[~,size_ltpks_cal]=size(ltpks_locs_cal);
[~,size_ltvls_cal]=size(ltvls_locs_cal);
[~,size_rtpks_cal]=size(rtpks_locs_cal);
[~,size_rtvls_cal]=size(rtvls_locs_cal);


if size_ltpks==size_ltpks_cal && size_ltvls==size_ltvls_cal && size_rtpks==size_rtpks_cal && size_rtvls==size_rtvls_cal % if the number of peaks and valleys are the same for calculation and experiemnt
    difference(i,dx)=sumsqr(ltpks_locs-ltpks_locs_cal)+sumsqr(ltvls_locs-ltvls_locs_cal)+sumsqr(rtpks_locs-rtpks_locs_cal)+sumsqr(rtvls_locs-rtvls_locs_cal);
end


 end
 
 
end
%% Surface plot of the difference function
figure
surf(diff_pks);

[min_val,idx]=min(difference(:));
[thickness,shift]=ind2sub(size(difference),idx);


figure
surf(difference);
title('See the minima for candidate thicknesses')
ylabel('thickness') 

disp('if the voltage is 200 kV, then the thickness may be (refer to JEMS simulation and the surface plot to verify this)')
disp(thickness)
disp('the improvement in the thickness measurement may rely on new alogrithms, which may invoke deep learning')











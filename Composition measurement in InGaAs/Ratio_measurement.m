%%The code is written by Dr. Yueming Guo to get the ratio measurement for his composition measurement in InGaAs 
%% Run this section first before proceeding
clear;
clc;
close all;
P_rotated=imread('centre from bottom 11.tif');
[totrow,totcol]=size(P_rotated);
Q_rotated=zeros(totrow+400,totcol+400);
Q_rotated(1:totrow,1:totcol)=P_rotated(1:totrow,1:totcol);
P_rotated=Q_rotated;
load('I400_9mrad.mat');
% Rough knowledge about the In composition, X_In, you can choose from X_In=0.00, 0.25, 0.50, 0.75, 1.00
%setting the scalebar
X_In=0.50;% roughly estimated composition
comp_number=3;% #1: 0   #2: 0.25   #3: 0.50  #4:  0.75  #5:  1.00 
scalebar_initial=36.9; %pixels/(1/nm) keep 3 digits at most @240mm Titan 36.9 @ 300mm F20 this is 70.9 @ Titan this is 37.1/37.2 @Titan Aug22 35.8
%%
Displayed_P_rotated=log10(1+0.0005*abs(P_rotated));
figure
imshow(Displayed_P_rotated);
%%
%find the maximum point from the pattern to help location of the central
%disc
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

%% Find the centre of the 000 disc and find the diametre of the disc
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
%% now draw a rectangle that covers the 400 kicuchi bands and record the top-left and bottom-right corners of the
%rectangle

x_br=round(col_centre+g200length+0.5*diametre);
y_br=round(row_centre-g400length);
x_tl=round(x_br-g400length-diametre);
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

%% Take the profile of the summed kikuchi band
Kikuchi_profile=sum(crop_kikuchi);
Kikuchi_profile=smooth(Kikuchi_profile);
figure
plot(1:x_max,Kikuchi_profile);
smoothed_Kikuchi=smooth(Kikuchi_profile);
%smoothed_Kikuchi=smoothed_Kikuchi';
[size_smoothed,~]=size(smoothed_Kikuchi);
figure
plot(1:size_smoothed,smoothed_Kikuchi);
%% Take the first-order derivative
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
%% Insert here

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



%%
Kikuchi_width=X_R_Kikuchi-X_L_Kikuchi;
centre_Kikuchi=0.5*(X_R_Kikuchi+X_L_Kikuchi);
X_L=centre_Kikuchi-Kikuchi_width/4;
X_R=centre_Kikuchi+Kikuchi_width/4;
X_L=X_L+x_tl-1;
X_R=X_R+x_tl-1;
centre_Kikuchi=centre_Kikuchi+x_tl-1;
%setting the scalebar
% scalebar_initial=41.2; %pixels/(1/nm)
scalebar=scalebar_initial*10; %pixels/(1/A)
E0=200*1e3; % accelrating voltage
lambda=12.2643/sqrt(E0*(1+0.978476*1e-6*E0)); % the relativiscally corrected wavlength
scalebar=scalebar*1e-3/lambda; %pixels/mrad
scalebar=round(scalebar); %pixels(integer)/mrad
pixel_1mrad=scalebar;
X_L_min=round(X_L-scalebar);
X_L_max=round(X_L+scalebar);
X_R_min=round(X_R-scalebar);
X_R_max=round(X_R+scalebar);

%% Get the line profile along the X_L and X_R
% row_min=300;
% row_max=700;

left_profile=P_rotated([row_min:row_max],[X_L_min:X_L_max]);
right_profile=P_rotated([row_min:row_max],[X_R_min:X_R_max]);

left_profile=sum(left_profile,2);
right_profile=sum(right_profile,2);
figure
plot(row_min:row_max,left_profile);
hold on
plot(row_min:row_max,right_profile);

%% Get the maximum of the profile

[max_num, ~]=max(left_profile);
[Y_L]=ind2sub(size(left_profile),find(left_profile==max_num));

[max_num, ~]=max(right_profile);
[Y_R]=ind2sub(size(right_profile),find(right_profile==max_num));

Y_L_min=Y_L-scalebar;
Y_L_max=Y_L+scalebar;

Y_R_min=Y_R-scalebar;
Y_R_max=Y_R+scalebar;

I_L=sum(left_profile(Y_L_min:Y_L_max));
I_R=sum(right_profile(Y_R_min:Y_R_max));

% calculate the background

Y_L_bg_min_upper=Y_L-5*scalebar;
Y_L_bg_max_upper=Y_L-3*scalebar;
I_L_bg_upper=sum(left_profile(Y_L_bg_min_upper:Y_L_bg_max_upper));

Y_L_bg_min_lower=Y_L+3*scalebar;
Y_L_bg_max_lower=Y_L+5*scalebar;
I_L_bg_lower=sum(left_profile(Y_L_bg_min_lower:Y_L_bg_max_lower));

Y_R_bg_min_upper=Y_R-5*scalebar;
Y_R_bg_max_upper=Y_R-3*scalebar;
I_R_bg_upper=sum(left_profile(Y_R_bg_min_upper:Y_R_bg_max_upper));


Y_R_bg_min_lower=Y_R+3*scalebar;
Y_R_bg_max_lower=Y_R+5*scalebar;
I_R_bg_lower=sum(left_profile(Y_R_bg_min_lower:Y_R_bg_max_lower));

if (I_L_bg_upper<0)||(I_L_bg_lower<=0)|| (I_R_bg_upper<=0)||(I_R_bg_lower<=0)
    disp('the data seems to be very noisy, try a CBED pattern recorded under a longer exposure');
else
    disp('the data is good(not noisy)');
end
% Calculate the ratio

ratio=(I_L-0.5*(I_L_bg_upper+I_L_bg_lower))/(I_R-0.5*(I_R_bg_upper+I_R_bg_lower));

if (X_In>0.165)&&(ratio<1)
    ratio=1/ratio;
end
if (X_In<0.165)&&(ratio>1)
    ratio=1/ratio;
end
disp('the ratio is');
disp(ratio);


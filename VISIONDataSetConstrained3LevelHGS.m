
% % Constrained 3 Level HGS method.
% This code was arranged for quickly running pre-calculated form of VISION dataset shared
%       under the same GitHub repository.
% It takes 3 arguments as input:
%     -visionPath: The path of the VISION folder in the same GitHub Repository. This folder
%	  contains pre-calculated 3-Level HGS transforms, corresponding video frames, and
%          camera fingerprints.
%
%      -isPreCalculated: If it is set 0, code will calculate transform for frames and fingerprint.
%          It can be usable for running code for other 2 datasets or new videos.
%          If it is not set 0, code will use pre-calculated transforms in the VISION folder.
%
%      -isWeightExist: Due to the size of the weighting masks, they were not provided in the
%	 repository. If this parameter set 0, code will skip weighting, otherwise, it will read
%  	 corresponding mat files for weighting.
% As an output it will write to the console number of correctly verified videos and PCE
%       threshold
% %
function [] = VISIONDataSetConstrained3LevelHGS(visionPath,isPreCalculated,isWeightExist)

% Videos in VISION folder have 3 parent folders, one of them for stabilization
%   level,one of them for their identifier and the last one for non-match
%   or match test (true positive (tp) and false positive (fp))

videos = dir(strcat(visionPath,'*/*/*p'));


for i=1:length(videos)
    tic
    videoPath=videos(i).folder;
    videoResPath = strcat(videos(i).folder, '/', videos(i).name);
    fprintf('Running video is %s\n', videoResPath)
    if isPreCalculated == 0
        clear Fp_cam_app
        load(strcat(videoResPath,'/fp.mat'))
        calculate_Transforms(videoPath,Fp_cam_app,videoResPath)
    end
    
    
    % all the constrained 3-levels HGS are made by using
    % n_sub = 2, pce_sub = 2 and pce_vld = 42 in
    % transform validation step.
    find_Result(videoPath,videoResPath,2,2,42,isWeightExist);
    toc
end

result(visionPath)

end

function [] = calculate_Transforms(videoPath,Fp_cam_app,videoResPath)

frames = dir(strcat(videoPath,'/*.jpg'));
for i = 1:min(10, length(frames))
    frames_path(i).name = strcat(frames(i).folder,'/',frames(i).name);
end


crop_size = 250;
ncc_range = 65;
y_size = (size(Fp_cam_app,1)+2)/2;
x_size = (size(Fp_cam_app,2)+2)/2;
transform = [];
transform_1 = [];
transform_2 = [];

for i=1:length(frames_path)
    
    
    frame = imread(frames_path(i).name);
    Noisex_s = NoiseExtractFromImage(frame,2);
    Noisex_s = WienerInDFT(Noisex_s,std2(Noisex_s));
    
    disp(i)
    
    Fp_camera_1 = Fp_cam_app((y_size-crop_size-ncc_range):(y_size-1+crop_size+ncc_range),(x_size-crop_size-ncc_range):(x_size-1+crop_size+ncc_range));
    Noisex_ref = Noisex_s((y_size-crop_size):(y_size-1+crop_size),(x_size-crop_size):(x_size-1+crop_size));
    
    transform=try_transform(Fp_camera_1,Noisex_ref,4,transform,i,ncc_range,zeros(1,11));
    toc;
    
    %% trying transforms of coarse sample points was finished
    
    
    % find the best 5 coarse sample points transform for the current frame
    temp_transforms4frame = transform(transform(:,1)==i,:);
    temp_transforms4frame = sortrows(temp_transforms4frame,2);
    best_sample_transforms = temp_transforms4frame((end-4:end),:);
    
    clear temp_transforms4frame
    
    for point_index=1:5
        
        transform_center = best_sample_transforms(point_index,:);
        transform_1=try_transform(Fp_camera_1,Noisex_ref,2,transform_1,i,ncc_range,transform_center);
    end
    
    
    %% trying transforms of finer sample points was finished
    
    % find the best 5 finer sample points transform for the current frame
    temp_transforms4frame = transform_1(transform_1(:,1)==i,:);
    temp_transforms4frame = sortrows(temp_transforms4frame,2);
    best_sample_transforms = temp_transforms4frame((end-4:end),:);
    
    clear temp_transforms4frame
    
    for point_index=1:5
        
        transform_center = best_sample_transforms(point_index,:);
        transform_2=try_transform(Fp_camera_1,Noisex_ref,1,transform_2,i,ncc_range,transform_center);
    end
    
    
    
    save(strcat(videoResPath,'/all_transforms.mat'),'transform','transform_1','transform_2');
end

disp('all transforms finished ...');
end

function [transform] = try_transform(Fp_camera_1,Noisex_ref,step,transform,frame_id,ncc_range,transform_center)
transform_count = size(transform,1);
transform_center_point=reshape(transform_center(4:11),[2 4]);
transform_center_point=transform_center_point';
ncc_range2=2*ncc_range;
for x1=-1:1
    for x2=-1:1
        for x3=-1:1
            for y1=-1:1
                for y2=-1:1
                    for y3=-1:1
                        
                        movingPoints = [1, 1; size(Fp_camera_1,2), 1; 1, size(Fp_camera_1,1); size(Fp_camera_1,2), size(Fp_camera_1,1)];
                        fixedPoints= [1+step*x1, 1+step*y1;...
                            size(Fp_camera_1,2)+step*x2, 1+step*y2;...
                            1+step*x3, size(Fp_camera_1,1)+step*y3;...
                            size(Fp_camera_1,2), size(Fp_camera_1,1)];
                        
                        fixedPoints =fixedPoints + transform_center_point;
                        
                        tform = fitgeotrans(movingPoints,fixedPoints,'projective');
                        
                        Noisex_t = imwarp(Noisex_ref,tform);
                        
                        
                        try
                            C = circxcorr2(Fp_camera_1,Noisex_t);
                            detection_1_fp = PCE2(C,[ncc_range2 ncc_range2]);
                            res = detection_1_fp.PCE;
                        catch
                            disp('PCE calculation error')
                            res = 0;
                        end
                        transform_count=transform_count+1;
                        transform(transform_count,1) = frame_id;
                        transform(transform_count,2) = res;
                        transform(transform_count,4) = step*x1;
                        transform(transform_count,5) = step*y1;
                        transform(transform_count,6) = step*x2;
                        transform(transform_count,7) = step*y2;
                        transform(transform_count,8) = step*x3;
                        transform(transform_count,9) = step*y3;
                        transform(transform_count,10) = 0;
                        transform(transform_count,11) = 0;
                        transform(transform_count,4:11)=transform(transform_count,4:11)+transform_center(4:11);
                    end
                end
            end
        end
    end
end
end

function [pceRes]=find_Result(videoPath,save_name,pce_sub,nsub,pce_glob,isWeightExist)

load(strcat(save_name,'/all_transforms.mat'));
load(strcat(save_name,'/fp.mat'));
[best20]=find_best_20(transform_2);

frames = dir(strcat(videoPath,'/*.jpg'));
%weights should have the same name with their frames except for their extensions.
for i = 1:min(10, length(frames))
    frames_path(i).name = strcat(frames(i).folder,'/',frames(i).name);
    frame_name=strsplit(frames(i).name,'.jpg');
    weigths_path(i).name = strcat(frames(i).folder,'/',frame_name{1},'.mat');
end


crop_size = 250;
ncc_range = 65;
y_size = (size(Fp_cam_app,1)+2)/2;
x_size = (size(Fp_cam_app,2)+2)/2;

limit = min(length(frames_path),max(best20(:,1)));


for i = 1:limit
    
    best20=createSubBlocks(best20,i,frames_path,pce_glob,Fp_cam_app);
    
end



Fp_camera_1 = Fp_cam_app((y_size-crop_size-ncc_range):(y_size-1+crop_size+ncc_range),(x_size-crop_size-ncc_range):(x_size-1+crop_size+ncc_range));
for i=1:5
    aggregateFP{i} = 0;
    aggregateFPD{i} = 1;
end

for i = 1:min(size(frames_path,2),size(best20,1))
    
    candidateTransformDetail = best20(((i-1)*20+1):(i*20),:);
    candidateTransformDetail(:,18) = sum(candidateTransformDetail(:,14:17)>pce_sub,2);
    validTransforms = candidateTransformDetail(candidateTransformDetail(:,18)>nsub,:);
    
    maxPceGlob=max(validTransforms);
    
    
    if size(validTransforms,1)>4 && maxPceGlob(2)>pce_glob
        frame = imread(frames_path(i).name);
        Noisex_s = NoiseExtractFromImage(frame,2);
        Noisex_s = WienerInDFT(Noisex_s,std2(Noisex_s));
        
        if isWeightExist == 1
            load(weigths_path(i).name)
        else
            weight = ones(size(Noisex_s));
        end
        Noisex_ref = Noisex_s((y_size-crop_size):(y_size-1+crop_size),(x_size-crop_size):(x_size-1+crop_size));
        weight_ref = weight((y_size-crop_size):(y_size-1+crop_size),(x_size-crop_size):(x_size-1+crop_size));
        
        for candidates=0:4
            validTransform=validTransforms((end-candidates),:);
            
            transform_point=reshape(validTransform(1,4:11),[2 4]);
            transform_point=transform_point';
            movingPoints = [1, 1; size(Fp_camera_1,1), 1; 1, size(Fp_camera_1,1); size(Fp_camera_1,2), size(Fp_camera_1,1)];
            fixedPoints= movingPoints + transform_point;
            
            tform = fitgeotrans(movingPoints,fixedPoints,'projective');
            [Noisex_t] = imwarp(Noisex_ref,tform);
            [Weight_t] = imwarp(weight_ref,tform);
            
            C = circxcorr2(Fp_camera_1,Noisex_t);
            detection_1_fp = PCE2(C,[130 130]);
            fw_best_crop = 132-detection_1_fp.PeakLocation;
            
            Noisex_crop = Noisex_t(fw_best_crop(1):fw_best_crop(1)+365-1,fw_best_crop(2):fw_best_crop(2)+365-1);
            Weight_crop = Weight_t(fw_best_crop(1):fw_best_crop(1)+365-1,fw_best_crop(2):fw_best_crop(2)+365-1);
            
            if isWeightExist ~= 1
                Weight_crop(Weight_crop==0)=1;
            end
            aggregateFP{candidates+1} = aggregateFP{candidates+1}+Noisex_crop.*Weight_crop;
            aggregateFPD{candidates+1} = aggregateFPD{candidates+1}+Weight_crop;
        end
    end
end


Fp_camera_1_t =Fp_camera_1(130:630, 130:630);

for candidates=1:5
    if aggregateFP{candidates} == 0
        pceRes(candidates) = 0;
    else
        C = circxcorr2(Fp_camera_1_t,aggregateFP{candidates}./aggregateFPD{candidates});
        detection_1_fp = PCE2(C,[3 3]);
        
        pceRes(candidates) = detection_1_fp.PCE;
    end
end

save(strcat(save_name,'/Top5Results.mat'),'pceRes');
end

function [best20] = createSubBlocks(best20,i,frames_path,pce_glob,Fp_cam_app)

crop_size = 250;
ncc_range = 65;
y_size = (size(Fp_cam_app,1)+2)/2;
x_size = (size(Fp_cam_app,2)+2)/2;

maxValue=max(best20((i-1)*20+1:i*20,2));
if maxValue < pce_glob
    best20((i-1)*20+1:i*20,13:17) = zeros(20,5);
    return
end
frame = imread(frames_path(i).name);
Noisex_s = NoiseExtractFromImage(frame,2);
Noisex_s = WienerInDFT(Noisex_s,std2(Noisex_s));

disp(i)

Fp_camera_1 = Fp_cam_app((y_size-crop_size-ncc_range):(y_size-1+crop_size+ncc_range),(x_size-crop_size-ncc_range):(x_size-1+crop_size+ncc_range));
Noisex_ref = Noisex_s((y_size-crop_size):(y_size-1+crop_size),(x_size-crop_size):(x_size-1+crop_size));

for x1=1:20
    transform_point=reshape(best20((i-1)*20+x1,4:11),[2 4]);
    transform_point=transform_point';
    movingPoints = [1, 1; size(Fp_camera_1,1), 1; 1, size(Fp_camera_1,1); size(Fp_camera_1,2), size(Fp_camera_1,1)];
    fixedPoints= movingPoints + transform_point;
    tform = fitgeotrans(movingPoints,fixedPoints,'projective');
    
    [Noisex_t] = imwarp(Noisex_ref,tform);
    
    C = circxcorr2(Fp_camera_1,Noisex_t);
    detection_1_fp = PCE2(C,[130 130]);
    
    best20((i-1)*20+x1,13) = detection_1_fp.PCE;
    
    fw_best_crop = 132-detection_1_fp.PeakLocation;
    
    Fp_camera_1_t =Fp_camera_1(130:630, 130:630);
    Noisex_crop = Noisex_t(fw_best_crop(1):fw_best_crop(1)+365-1,fw_best_crop(2):fw_best_crop(2)+365-1);
    
    C = circxcorr2(Noisex_crop(1:crop_size,1:crop_size),Fp_camera_1_t(1:crop_size,1:crop_size));
    detection_1_fp = PCE2(C,[2 2]);
    best20((i-1)*20+x1,14) = detection_1_fp.PCE;
    C = circxcorr2(Noisex_t(1:crop_size,end-crop_size:end),Fp_camera_1_t(1:crop_size,end-crop_size:end));
    detection_1_fp = PCE2(C,[2 2]);
    best20((i-1)*20+x1,15) = detection_1_fp.PCE;
    C = circxcorr2(Noisex_t(end-crop_size:end,1:crop_size),Fp_camera_1_t(end-crop_size:end,1:crop_size));
    detection_1_fp = PCE2(C,[2 2]);
    best20((i-1)*20+x1,16) = detection_1_fp.PCE;
    C = crosscorr(Noisex_t(end-crop_size:end,end-crop_size:end),Fp_camera_1_t(end-crop_size:end,end-crop_size:end));
    detection_1_fp = PCE2(C,[2 2]);
    best20((i-1)*20+x1,17) = detection_1_fp.PCE;
    
end
end
function [best20]=find_best_20(transform_2)
sizes = max(transform_2(:,1));

for i=1:sizes
    max_indis = find(transform_2(:,1)==i);
    sonucc = transform_2(max_indis,:);
    sonucc = sortrows(sonucc,2);
    best20((i*20-19):(i*20),:) = sonucc((end-19):end,:);
    clear sonucc
end
end

function [] = result(visionPath)
trueFiles = dir(strcat(visionPath,'*/*/tp/Top5Results.mat'));
falseFiles = dir(strcat(visionPath,'*/*/fp/Top5Results.mat'));

trues=loadFiles(trueFiles);
falses=loadFiles(falseFiles);
pce_threshold=max(falses(:,3));
fprintf('The number of correctly attributed videos: %d\n',sum(trues(:,3)>pce_threshold));
fprintf('PCE threshold: %d\n',pce_threshold)

end

function [results] =loadFiles(files)
for i=1:length(files)
    clear pceRes
    load(strcat(files(i).folder,'/',files(i).name))
    results(i,:)=sort(pceRes);
end
end


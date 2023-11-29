%Doing yeast cell segmentation from 2-channel images (gfp and rfp)
infld='ImagingData/';
outfld = 'ImagingData/Masks/';
flist=dir(infld);
mkdir('ImagingData/Masks/');
flist=flist(3:end);
for s=1:length(flist)
     s
    if(findstr(flist(s).name,'ome.tiff'))
        fnme=fullfile(infld, flist(s).name)
        gfpim = double(imread(fnme,1));
        rfpim = double(imread(fnme,2));
        matfnme=[outfld,flist(s).name(1:(end-9)),'.mask.mat'];
        %gfpim = double(imread(fnme,1));
        %rfpim = double(imread(fnme,2));
        %if(sum(gfpim(:)) ~=0) 
        cellbndry = gfpim+rfpim;
        cellbndry_tmp=imread(fnme,1)+imread(fnme,2);
        %figure;subplot(1,2,1); imshow(imadjust(cellbndry_tmp));
        [f,xi] = ksdensity(cellbndry(:));
        df=diff(f);
        [pks,locs] = findpeaks(f);
        flg=0;cnt=1;
        for j=1:size(df,2)%locs(1,1):size(df,2)
            if(flg==0)thresh=xi(j)+1;flg=1;else thresh=xi(j); end
            celledge = imbinarize(cellbndry/max(cellbndry(:)),thresh/max(cellbndry(:)));
            celledge1=bwfill(celledge,'holes');
            D = -bwdist(~celledge1);%figure;imshow(D,[]);
            %Getting markers
            markers = imextendedmin(D,2);%figure;imshow(markers,[]);
            %Imposing minima at markers pixels
            D2 = imimposemin(D,markers);
            %Performing watershed transformation
            Ld2 = watershed(D2);
            Mask1 = celledge1;Mask1(Ld2 == 0) = 0;
            if(cnt == 1)Mask_tmp=Mask1;else Mask_tmp(:,:,cnt)=Mask1;end
            cnt=cnt+1;
            %figure;imshow(Mask1);
            if (df(1,j) >=0 && cnt >15)
            	break
            end
        end
        Mask_combine=zeros(size(Mask1));
        markers_tot=zeros(size(Mask1));
        rtot=[];ctot=[];figure;
        for i=1:(cnt-1)
            subplot(1,2,2);
            if(i<(cnt-1))
                Ltmp1=bwlabel(Mask_tmp(:,:,i+1));            
                imshow(label2rgb(Ltmp1,'jet',[1 1 1],'shuffle'),[]);
            else Ltmp1 = zeros(size(Mask1)); 
                imshow(Ltmp1,[]);
            end
            Ltmp=bwlabel(Mask_tmp(:,:,i));
            subplot(1,2,1);imshow(label2rgb(Ltmp,'jet',[1 1 1],'shuffle'),[]);
            title([num2str(i),'/',num2str((cnt-1))]);
            for j = 1:size(ctot,1)
                subplot(1,2,1);hold on;plot(ctot(j,1),rtot(j,1),'*k');
            end
            [c r]=getpts();
            for j = 1:size(c,1)
                labelval=Ltmp(round(r(j)),round(c(j)));
                Mask_combine(Ltmp==labelval)=1;
                markers_tot(round(r(j)),round(c(j)))=1;
            end
            ctot=[ctot;c];rtot=[rtot;r];
        end
        hold off;
        %Start resegmenting final mask
        %se=strel('disk',2);a1=imclose(Mask_combine,se);
        D = -bwdist(~Mask_combine);D2 = imimposemin(D,markers_tot);
        Ld2 = watershed(D2);Mask_combine(Ld2 == 0) = 0;
        %title('Outlined Original Image')
        Ltmp=bwlabel(Mask_combine);
        overlay1 = imoverlay(imadjust(cellbndry_tmp), BWoutline, [.3 1 .3]);
        subplot(1,2,1);imshow(overlay1);
        subplot(1,2,2);imshow(label2rgb(Ltmp,'jet',[1 1 1],'shuffle'),[]);
        %Removing cells on border
        %Mask = imclearborder(Mask_combine);
        L=bwlabel(Mask_combine);
        %Removing the stray pixels.
         LabelProp = regionprops(L,'Area');
         Areas=[LabelProp.Area];
         for i=1:size(Areas,2)
             if(Areas(1,i) < 500)
                 L(L==i)=0;
             end
         end
         %figure;imshow(label2rgb(L),[]);
         Mask_final=L;
         Mask_final(L>0)=1;
         save(matfnme,'Mask_final');
   end
   clearvars fnme matfnme;
end

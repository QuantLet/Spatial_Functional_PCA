[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **SFPCA_divide_data** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

﻿
Name of Quantlet: SFPCA_divide_data

Published in: Spatial_Functional_PCA

Description: 'Divides the 91*109*91 voxels data into 30*35*30 (or any other dimensions) cubes for each subject, related to the paper Spatial Functional Principal Component Analysis with Applications to Brain Image Data by Yingxing Li, Chen Huang and Wolfgang Härdle.'

Keywords : 
- fmri
- pca
- smoothing
- basis
- factor analysis

Author:          Chen Huang

Submitted:        Wednesday, June 20 2018 by Chen Huang

```

### MATLAB Code
```matlab

  subjects = [1,2,3,4,5,6,8,9,10,11,12,15,16,17,18,19,21];
  tic;
  load('scanidx')
  for s=1:length(subjects)
    subject = subjects(s);
    filename = ['fmridata',num2str(subject),'cnew'];
    load(filename)
%     filename = ['scanidx',num2str(subject)];
%     load(filename)
    for q=1:81
      temp_q = temp(:,:,:,scan_idx(q,:)); 
      filename = ['voxels_',num2str(subject),'_',num2str(q)];
      save(filename, 'temp_q')
    end
  end
  toc;
```

automatically created on 2018-09-04
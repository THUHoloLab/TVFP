## Total variation regularization-based Fourier ptychography (TVFP)<br><sub>Matlab implementation of the TVFP / eTVFP</sub>

<p align="center">
<img src="https://github.com/THUHoloLab/TVFP/blob/master/figures/fig1.png" width="500">
</p>
  
**Resolution enhancement of long-range imaging with sparse apertures**<br>
Jiachen Wu, Feng Yang, Liangcai Cao*<br>
- [Paper link](https://www.sciencedirect.com/science/article/abs/pii/S0143816622001208) for TVFP / eTVFP <a href="https://www.sciencedirect.com/science/article/abs/pii/S0143816622001208"></a>.
- [Paper download](https://github.com/THUHoloLab/TVFP/raw/master/figures/OLEN-jiachen.pdf) for TVFP / eTVFP <a href="https://github.com/THUHoloLab/TVFP/raw/master/figures/OLEN-jiachen.pdf"></a>.
- Contact: wjc18@mails.tsinghua.edu.cn
- Contact: clc@tsinghua.edu.cn

## Release notes
- demo_TVFP_resChart.m  - demonstrate the TVFP algorithm to reconstruct amplitude-only objects (resolution chart)  
- demo_TVFP_complex.m   - demonstrate the TVFP algorithm to reconstruct complex objects  
- demo_eTVFP.m          - demonstrate the eTVFP algorithm which could estimate pupil function  
- demo_realData.m       - demonstrate the reconstruction for real data  

## Performance evaluation
Reconstructions for different overlapping ratio. (a) Reconstructed images under different overlapping ratio with Gaussian noise using different algorithms (Gauss-Newton, TPWFP, FPSR, and TVFP). (b) PSNR for reconstructed images from measurements without noise. (c) PSNR for reconstructed images from measurements with 20 dB SNR Gaussian noise.

<p align="center">
<img src="https://github.com/THUHoloLab/TVFP/blob/master/figures/fig2.png" width="700">
</p>

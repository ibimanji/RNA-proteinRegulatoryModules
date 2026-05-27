conda环境配置：

conda activate pystar310
conda install -y -c conda-forge numpy scipy matplotlib networkx scikit-image scikit-learn tifffile zarr pyjnius blessed
pip install ashlar
conda install pandas numpy scikit-image
conda install -c conda-forge "zarr<3.0.0" -y
conda install -c conda-forge "numcodecs<0.14" -y




sbatch 【run_stitch.sh】 \
    【输入文件夹，其中每个Position包含cell_center.csv, remain_reads_raw.csv】 \
    【输出文件夹】\
    【MAF文件路径】 \
    【xy pixelsize - μm】 \
    【DAPI通道】 \
    【FOV旋转角度 - 顺时针】 \
    【PI/DAPI图像文件夹，可选，例如02_registration/IF/PI】 \
    【MAF起始PositionID，可选】 \
    【MAF结束PositionID，可选】 \
    【点/细胞坐标旋转角度 - 顺时针，可选，默认等于FOV旋转角度】 \
    【局部overlay QC数量，可选，默认8】

示例：
sbatch run_stitch.sh \
    /gpfs/share/home/2301920002/labShare/2301920002/ADdecon/14mWT_14mFAD/14mWT/03_segmentation \
    /gpfs/share/home/2301920002/labShare/2301920002/ADdecon/14mWT_14mFAD/14mWT/04_stitch \
    /gpfs/share/home/2301920002/labShare/2301920002/ADdecon/14mWT_14mFAD/output_2025-11-01_20251101A2B2_14mWT_14mAD.maf \
    0.0946 \
    ch00 \
    90 \
    /gpfs/share/home/2301920002/labShare/2301920002/ADdecon/14mWT_14mFAD/14mWT/02_registration/IF/PI \
    1 \
    127 \
    0 \
    8

输出包含：
    FOV 对齐后的物理坐标 (μm)：registered_tilecoord.csv
    细胞中心拼接坐标：cell_centerouter.csv
    转录本拼接坐标：remain_readsouter.csv
    可视化图像：cell_reads_profile.png
    各通道拼接结果：chxx_stitched.ome.tif/stitched.ome.tif

备注：
    如果提供MAF起始/结束PositionID，本地Position文件夹会从Position001开始重新匹配：
        MAF起始Position -> 本地Position001，例如MAF 1-127匹配本地1-127，MAF 128-254也匹配本地1-127。
    MAF起始/结束PositionID必须是连续范围。
    如果使用02_registration/IF/PI中的原始PI图像，且原始PI需要顺时针90度后才能和MAF对应，图像旋转用90。
    clustermap的remain_reads_raw.csv已经和顺时针90度后的PI坐标对应；这种情况下点/细胞坐标旋转用0。
    程序会输出local_overlay_qc/PositionXXX_overlay.png，用来检查单个FOV里点/细胞中心是否和旋转后的PI图像对齐。
    如果某个Position没有remain_reads_raw.csv，会用空表代替，不中断图像拼接。

import os
import numpy as np
import skimage.io
from ashlar import reg
from ashlar.reg import EdgeAligner

class DirectTiffReader(reg.PlateReader):
    def __init__(self, pos_coords, channel_str, pixel_size, rotation=0):
        self.pos_coords = pos_coords
        self.channel_str = channel_str
        self.real_pixel_size = pixel_size
        self.k_val = -(rotation // 90)
        
        self.files = []
        for entry in self.pos_coords:
            pos_path = entry[0]
            image_path = entry[2] if len(entry) > 2 else None
            if image_path:
                self.files.append(image_path)
                continue
            if self.channel_str == "":
                tifs = [f for f in os.listdir(pos_path) if f.endswith('.tif')]
            else:
                tifs = [f for f in os.listdir(pos_path) if self.channel_str in f and f.endswith('.tif')]
                
            if not tifs:
                raise FileNotFoundError(f"Cannot find valid TIF for channel '{self.channel_str}' in {pos_path}")
            self.files.append(os.path.join(pos_path, tifs[0]))
        
        img = skimage.io.imread(self.files[0])
        if img.ndim == 3: img = np.max(img, axis=0)
        img_rot = np.rot90(img, k=self.k_val)
        
        self.tile_size_val = np.array(img_rot.shape)
        self.pixel_dtype = img.dtype
        
    @property
    def metadata(self):
        class Meta: pass
        m = Meta()
        m.num_images = len(self.files)
        m.num_channels = 1
        m.pixel_size = self.real_pixel_size
        
        m.positions = np.array([
            [entry[1][1] / self.real_pixel_size, entry[1][0] / self.real_pixel_size]
            for entry in self.pos_coords
        ])
        
        m.size = self.tile_size_val
        m.filename = self.files
        m.origin = np.array([0, 0])
        m.pixel_dtype = self.pixel_dtype
        return m
        
    def read(self, series, c):
        img = skimage.io.imread(self.files[series])
        if img.ndim == 3: img = np.max(img, axis=0)
        return np.rot90(img, k=self.k_val)

def align_reference(pos_coords, ref_channel, pixel_size, rotation=0, max_shift_um=150.0):
    reader = DirectTiffReader(pos_coords, ref_channel, pixel_size, rotation)
    
    aligner = EdgeAligner(
        reader, 
        channel=0, 
        filter_sigma=1.0, 
        verbose=True,
        max_shift=max_shift_um
    )
    aligner.run()
    
    return aligner

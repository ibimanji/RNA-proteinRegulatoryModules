import os
from ashlar.reg import Mosaic, PyramidWriter
from align import DirectTiffReader

def stitch_all_channels(aligner, pos_coords, channels, output_dir, pixel_size, rotation=0):
    os.makedirs(output_dir, exist_ok=True)
    
    for ch in channels:
        channel_label = ch if ch else "Single Channel"
        print(f"\n[{channel_label}] Applying coordinates and stitching...")
        
        reader_ch = DirectTiffReader(pos_coords, ch, pixel_size, rotation)
        aligner.reader = reader_ch
        aligner._cache = {} 
        
        mosaic = Mosaic(aligner=aligner, shape=aligner.mosaic_shape, verbose=False)
        
        prefix = f"{ch}_" if ch else ""
        out_path = os.path.join(output_dir, f"{prefix}stitched.ome.tif")
        
        writer = PyramidWriter(
            mosaics=[mosaic], 
            path=out_path, 
            scale=2, 
            tile_size=1024, 
            verbose=False
        )
        writer.run()
        print(f"    Saved -> {out_path}")
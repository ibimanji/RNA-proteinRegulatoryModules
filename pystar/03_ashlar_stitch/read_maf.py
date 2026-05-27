import xml.etree.ElementTree as ET

def get_theoretical_coords(maf_path):
    # Parse MAF file and return physical coordinates mapped by PositionID
    tree = ET.parse(maf_path)
    points = tree.getroot().findall('.//XYZStagePointDefinition')
    
    maf_data = {}
    for pt in points:
        pos_id = int(pt.get('PositionID'))
        # Leica MAF coordinates are in meters, multiply by 1e6 to convert to um
        x = float(pt.get('StageXPos')) * 1e6
        y = float(pt.get('StageYPos')) * 1e6
        maf_data[pos_id] = {'x': x, 'y': y}
        
    if not maf_data:
        return {}
        
    # Find minimum values to zero-center the coordinates
    min_x = min(d['x'] for d in maf_data.values())
    min_y = min(d['y'] for d in maf_data.values())
    
    # Return a dictionary {pos_id: (x_um, y_um)}
    coords_dict = {
        pos_id: (d['x'] - min_x, d['y'] - min_y)
        for pos_id, d in maf_data.items()
    }
    
    return coords_dict
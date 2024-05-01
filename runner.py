import subprocess
import os
import xatlas
import trimesh
import numpy as np
from perlin_numpy import generate_perlin_noise_3d


# yoinked from https://gist.github.com/aminzabardast/cdddae35c367c611b6fd5efd5d63a326
'''
Write a Numpy array to a PFM file.
'''
def write_pfm(file, image, scale = 1):
  file = open(file, 'wb')

  color = None

  if image.dtype.name != 'float32':
    raise Exception('Image dtype must be float32.')

  if len(image.shape) == 3 and image.shape[2] == 3: # color image
    color = True
  elif len(image.shape) == 2 or len(image.shape) == 3 and image.shape[2] == 1: # greyscale
    color = False
  else:
    raise Exception('Image must have H x W x 3, H x W x 1 or H x W dimensions.')

  file.write(b'PF\n' if color else b'Pf\n')
  file.write(b'%d %d\n' % (image.shape[1], image.shape[0]))

  endian = image.dtype.byteorder

  if endian == '<' or endian == '=' and sys.byteorder == 'little':
    scale = -scale

  file.write(b'%f\n' % scale)

  image.tofile(file) 

def create_scene_file(input_file_path, output_file, transform, rotation=None, noise_texture_file=None):
    with open(output_file, 'w') as f_out:
        f_out.write(f'm {input_file_path}\n')
        f_out.write(f't {transform}\n')
        if rotation is not None:
            f_out.write(f'r {rotation}\n')
        if noise_texture_file is not None:
            f_out.write(f'n {noise_texture_file}')

def process_file(input_file_path, filename, rotation=None, noise_texture_slice=None):
    path_without_ext = input_file_path[:len(input_file_path) - 4]
    temp_file_path = path_without_ext + ".scene"
    transform = "0.0 0.0 0.0"
    if noise_texture_slice is None:
        noise_texture_file = None
    else:
        pass
    create_scene_file(input_file_path, temp_file_path, transform, rotation, noise_texture_file)
    
    output_png = path_without_ext + ".png"
    command = [
        "./pathtracer", 
        "-t", "72", 
        "-j", temp_file_path, 
        "-e", "../exr/outdoor_umbrellas_16k.exr", 
        "-f", output_png, 
        "-v", "1", 
        "-s", "1024", 
        "-m", "32",
        "-z", "1",
        "-n", "1",
        "-r", "1280", "720",
        "../dae/simple/empty.dae"]

    print(' '.join(command))
    
    subprocess.run(command, check=True)
    
def main():
    folder_path = input("specify working directory: ") ## todo 
    noise = input("noise? [y/n]: ") == 'y'
    spin = input("spin? [y/n]: ") == 'y'

    noise_texture = None
    rotation_rate = 0.0
    
    if noise:
        num_frames = sum([1 if s.endswith(".obj") for s in os.listdir(folder_path)])
        noise_texture = generate_perlin_noise_3d((1024, 1024, num_frames), (8, 8, 8))
    if spin:
        rotation_rate = float(input("how fast? degrees/s: ")) / 30.0

    rotation = 0.0
    # Iterate over all files in the folder
    i = 0
    for filename in os.listdir(folder_path):
        if filename.endswith(".obj"):  ## onyl want obj
            file_path = os.path.join(folder_path, filename)
            noise_texture_slice = noise_texture[:,:,i]
            process_file(file_path, filename, rotation=rotation)
            rotation += rotation_rate
            i += 1
            

if __name__ == "__main__":
    main()
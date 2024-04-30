import subprocess
import os

def create_scene_file(input_file_path, output_file, transform):
    with open(output_file, 'w') as f_out:
        f_out.write(f'{input_file_path}\n')
        f_out.write(f'{transform}')

def process_file(input_file_path, filename):
    path_without_ext = input_file_path[:len(input_file_path) - 4]
    temp_file_path = path_without_ext + ".scene"
    transform = "t 0.0 0.0 0.0" ## todo
    create_scene_file(input_file_path, temp_file_path, transform)
    
    output_png = path_without_ext + ".png"
    command = [
        "./pathtracer", 
        "-t", "8", 
        "-j", temp_file_path, 
        "-e", "../exr/little_paris_under_tower_2k.exr", 
        "-f", output_png, 
        "-v", "1", 
        "-s", "256", 
        "-m", "32"
        "-z", "1",
        "-r", "1280", "720",
        ""] ## todo
    
    subprocess.run(command, check=True)
    os.remove(temp_file_path)

def main():
    folder_path = input("specify working directory: ") ## todo 
    
    # Iterate over all files in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith(".obj"):  ## onyl want obj
            file_path = os.path.join(folder_path, filename)
            process_file(file_path, filename)

if __name__ == "__main__":
    main()
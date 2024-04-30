import subprocess
import os

def add_extra_line(input_file_path, output_file, transform):
    with open(input_file_path, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            f_out.write(line)
        f_out.write(transform)

def process_file(input_file_path, filename):
    temp_file_path = input_file_path[:len(input_file_path) - 4] + ".scene"
    transform = "specify transform here" ## todo
    add_extra_line(input_file_path, temp_file_path, transform)
    
    command = ["./pathtracer", "-t", "8", temp_file_path] ## todo
    
    subprocess.run(command, check=True)
    os.remove(temp_file_path)

def main():
    folder_path = "specify working directory" ## todo 
    
    # Iterate over all files in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith(".obj"):  ## onyl want obj
            file_path = os.path.join(folder_path, filename)
            process_file(file_path, filename)

if __name__ == "__main__":
    main()
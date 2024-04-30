import subprocess
import os

def add_extra_line(input_file_path, output_file, transform):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            f_out.write(line)
        # Add your extra line here
        f_out.write(transform)

def process_file(input_file_path, file_name):
    temp_file_path = input_file_path[:len(input_file_path) - 4] + ".scene"
    transform = "specify transform here"
    add_extra_line(input_file_path, temp_file_path, transform)
    
    command = ["args", "ags", "args", temp_file_path]
    
    subprocess.run(command, check=True)
    os.remove(temp_file)

def main():
    folder_path = "specify working directory"
    
    # Iterate over all files in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith(".obj"):  ## onyl want obj
            file_path = os.path.join(folder_path, filename)
            process_file(file_path, file_name)

if __name__ == "__main__":
    main()
import os

def read_txt_file(filename):
    with open(filename, 'r') as file:
        content = file.read()
    return content

def read_txt_to_list(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        lines = [line.strip() for line in lines]
    return lines

def append_line(filename, line):
    with open(filename, 'a') as file:
        file.write(line + '\n')

def copy_file(source, destination):
    with open(source, 'rb') as f_source:
        with open(destination, 'wb') as f_destination:
            # Read and write the file in chunks to conserve memory
            for chunk in iter(lambda: f_source.read(4096), b''):
                f_destination.write(chunk)
          
                
def create_task(site_id='line',
            run_path =r'C:/Users/xipeng/OneDrive - UGent/Desktop/create/run/',
            input_path = r'/scratch/gent/vo/000/gvo00074/vsc44253/phd/paper2/HPC/Response_slope/create/Inputs/',
            output_path = r'/scratch/gent/vo/000/gvo00074/vsc44253/phd/paper2/HPC/Response_slope/Outputs/'
            ):

    if not os.path.exists(output_path):
        os.mkdir(output_path)
        
        
    task_dir = run_path + site_id +'/'
    if not os.path.exists(task_dir):
        os.mkdir(task_dir)
    
        copy_file('./job_template.sh',task_dir+'./job.sh')
        
        with open(task_dir+'./job.sh', 'r+') as file:
            lines = file.readlines()
            
            lines.insert(1, '#PBS -N ' + site_id + '\n')
            file.seek(0)
            
            file.writelines(lines)

        append_line(task_dir+'job.sh', 'Rscript ' + task_dir +"run.R")


        
        copy_file('./run.R',task_dir+'run.R')

        with open(task_dir+'run.R', 'r+') as file:
            lines = file.readlines()
            
            lines.insert(0, 'index =' + '"' +site_id + '"\n')
            lines.insert(1, 'Input_path =' + '"' +input_path + '"\n')
            lines.insert(2, 'Output_path =' + '"' +output_path + '"\n')

            file.seek(0)
            
            file.writelines(lines)


if __name__ == '__main__':
    
    current_path = os.getcwd()

    run_path =  current_path +'/..' +r'/run/'
    input_path = current_path + r'/Inputs/'
    output_path = current_path +'/..' + r'/Outputs/'
    
    if not os.path.exists(run_path):
        os.mkdir(run_path)
        
    if not os.path.exists(run_path + 'run.sh'):
        copy_file('./run.sh',run_path + '/run.sh')

    site_info_file = 'run_list.txt'
    file_lines = read_txt_to_list(site_info_file)
    # file_lines = ["LA","LD","LKC","LNC",
    #              "B.thick","Conduit.d.","Crown.dia",#"Crown.H",
    #              "LPC","L.thick","Vcmax","R.depth","SLA",
    #              "SD","Stom.c","Max.H","WD"]
    for line in file_lines:
        print(line)
        create_task(
            site_id=line,
            run_path =run_path,
            input_path = input_path,
            output_path = output_path)
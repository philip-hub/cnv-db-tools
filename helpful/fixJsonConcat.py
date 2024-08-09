input_file_path = 'SJBALL013380_D1 3.json'
output_file_path = 'test_fix2.json'

with open(input_file_path, 'r') as file:
    lines = file.readlines()


new_lines = [line.rstrip('\n') + ',' + '\n' for line in lines]
new_lines[0] = '[' + new_lines[0]
new_lines[-1] = new_lines[-1].rstrip(',\n') + ']\n'

with open(output_file_path, 'w') as file:
    file.writelines(new_lines) #o(2n) of time 

print(f"Modified file saved as {output_file_path}")

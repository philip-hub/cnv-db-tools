def count_characters_in_file(file_path):
    try:
        with open(file_path, 'r') as file:
            content = file.read()
            character_count = len(content)
            return character_count
    except FileNotFoundError:
        print(f"The file {file_path} does not exist.")
        return None

# Replace 'your_file.txt' with the path to your file
file_path = 'master/master.tsv'
character_count = count_characters_in_file(file_path)

if character_count is not None:
    print(f"The file {file_path} contains {character_count} characters.")

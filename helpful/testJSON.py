import json

def is_valid_json(file_path):
    try:
        with open(file_path, 'r') as file:
            json.load(file)
        return True
    except ValueError as e:
        print(f"not valid JSON: {e}")
        return False
    except FileNotFoundError:
        print("file was not found.")
        return False
    except Exception as e:
        print(f"another error error occurred: {e}")
        return False

def fix_json(json_str):
    try:
        # Attempt to load the JSON string to check if it's valid
        return json.loads(json_str)
    except json.JSONDecodeError as e:
        print(f"JSONDecodeError: {e}")
        # If there's an error, let's try to fix some common issues
        json_str = json_str.replace("\'", "\"")  # Replace single quotes with double quotes
        json_str = json_str.replace("None", "null")  # Replace Python None with JSON null
        json_str = json_str.replace("True", "true")  # Replace Python True with JSON true
        json_str = json_str.replace("False", "false")  # Replace Python False with JSON false

        try:
            return json.loads(json_str)
        except json.JSONDecodeError as e:
            print(f"Failed to fix JSON: {e}")
            return None


file_path = 'test_fix2.json'

if is_valid_json(file_path):
    print("JSON is valid.")
else:
    print("JSON file is invalid.")




# with open(file_path, 'r') as file:
#     json_str = str(file.read())

# fixed_json = fix_json(json_str)

# if fixed_json is not None:
#     print("Fixed JSON:", fixed_json)
# else:
#     print("Could not fix JSON.")
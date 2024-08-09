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


file_path = 'SJBALL013380_D1 3.json'

if is_valid_json(file_path):
    print("JSON is valid.")
else:
    print("JSON file is invalid.")

import re

# Function to remove brackets from specified fields
def remove_brackets(text):
    patterns = [
        r'"lcvq":"\[',  # Matches "lcvq":"[
        r'"\]',         # Matches "]
        r'"v":"\[',     # Matches "v":"[
        r'"posv":"\[',  # Matches "posv":"[
        r'"ai'          # Matches "ai
    ]
    
    for pattern in patterns:
        text = re.sub(pattern, lambda x: x.group().replace('[', '').replace(']', ''), text)
        
    return text

# Read the JSON file as a text file
with open('your_json_file.json', 'r') as file:
    json_text = file.read()

# Remove brackets from specified fields
fixed_json_text = remove_brackets(json_text)

# Save the modified text to fixed_format.json
with open('fixed_format.json', 'w') as file:
    file.write(fixed_json_text)

print("Brackets removed and saved to fixed_format.json")

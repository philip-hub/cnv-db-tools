import numpy as np
from transformers import BertTokenizer
from PIL import Image, ImageOps
import os

def load_image(image_path, target_size=(256, 256)):
    # Load an image and resize it to the target size while maintaining the aspect ratio
    image = Image.open(image_path).convert("RGB")
    image.thumbnail(target_size, Image.ANTIALIAS)
    
    # Create a new image with white background
    new_image = Image.new("RGB", target_size, (255, 255, 255))
    
    # Paste the resized image onto the new image
    new_image.paste(image, ((target_size[0] - image.size[0]) // 2, (target_size[1] - image.size[1]) // 2))
    
    return np.array(new_image)

def tokenize_text(text, tokenizer, max_length=128):
    # Tokenize the text using the specified tokenizer
    tokens = tokenizer.encode(text, add_special_tokens=True, max_length=max_length, padding='max_length', truncation=True)
    return np.array(tokens)

def tokenize_image(image_array):
    # Flatten the image array to create a sequence of pixel values
    return image_array.flatten()

def combine_tokens(text_tokens, image_tokens=None):
    # Combine text tokens and image tokens
    if image_tokens is not None:
        return np.concatenate([text_tokens, image_tokens])
    return text_tokens

def read_text_from_file(file_path):
    # Read text from a file
    with open(file_path, 'r') as file:
        text = file.read().strip()
    return text

def save_tokens_to_npz(file_path, text_tokens, image_tokens=None):
    # Save the tokens to a .npz file
    if image_tokens is not None:
        np.savez_compressed(file_path, text_tokens=text_tokens, image_tokens=image_tokens)
    else:
        np.savez_compressed(file_path, text_tokens=text_tokens)

# Example usage
def tokenize_input(text_file_path, image_path=None, target_image_size=(256, 256), max_text_length=128, output_file_path='combined_tokens.npz'):
    # Initialize the tokenizer
    tokenizer = BertTokenizer.from_pretrained('bert-base-uncased')
    
    # Read and tokenize the text from file
    text = read_text_from_file(text_file_path)
    text_tokens = tokenize_text(text, tokenizer, max_length=max_text_length)
    print(f"Text tokens shape: {text_tokens.shape}")
    
    # If an image path is provided, load and tokenize the image
    if image_path:
        print(f"Loading image from: {image_path}")
        image_array = load_image(image_path, target_size=target_image_size)
        print(f"Image shape: {image_array.shape}")
        image_tokens = tokenize_image(image_array)
        print(f"Image tokens shape: {image_tokens.shape}")
    else:
        image_tokens = None
    
    # Combine the text tokens and image tokens
    combined_tokens = combine_tokens(text_tokens, image_tokens)
    
    # Save the combined tokens to an .npz file
    save_tokens_to_npz(output_file_path, text_tokens, image_tokens)
    
    return combined_tokens

# Directory paths
text_dir = 'model_info'  # Directory containing text files
image_dir = 'model_images'  # Directory containing image files
output_dir = 'encodings'  # Directory to save the .npz files

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Get list of text files
text_files = [f for f in os.listdir(text_dir) if f.endswith('.txt')]

# Supported image formats
image_formats = ['.png', '.jpg', '.jpeg']

# Iterate over text files and corresponding images
for text_file in text_files:
    text_file_path = os.path.join(text_dir, text_file)
    
    # Check for image in different formats
    image_file_path = None
    for fmt in image_formats:
        image_file_candidate = text_file.replace('info', 'image').replace('.txt', fmt)
        if os.path.exists(os.path.join(image_dir, image_file_candidate)):
            image_file_path = os.path.join(image_dir, image_file_candidate)
            break
    
    output_file_path = os.path.join(output_dir, text_file.replace('.txt', '.npz'))
    
    if image_file_path:
        print(f"Processing text file: {text_file_path} and image file: {image_file_path}")
    else:
        print(f"Processing text file: {text_file_path} without corresponding image")
    
    tokens = tokenize_input(text_file_path, image_path=image_file_path, output_file_path=output_file_path)
    print(f"Processed and saved {output_file_path}")

print("All files processed and saved.")

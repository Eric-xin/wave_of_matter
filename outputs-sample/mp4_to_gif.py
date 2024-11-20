import subprocess
import os

def mp4_to_gif(mp4_path, gif_path):
    # Use ffmpeg to convert mp4 to gif
    subprocess.run([
        'ffmpeg', 
        '-i', mp4_path, 
        '-vf', 'fps=10,scale=320:-1:flags=lanczos', 
        '-c:v', 'gif', 
        '-dpi', '100',
        gif_path
    ])

# get a list of all the mp4 files in the outputs directory
mp4_files = [f for f in os.listdir('outputs') if f.endswith('.mp4')]
# convert each mp4 file to a gif
for mp4_file in mp4_files:
    mp4_path = os.path.join('outputs', mp4_file)
    gif_path = os.path.join('outputs', mp4_file.replace('.mp4', '.gif'))
    mp4_to_gif(mp4_path, gif_path)
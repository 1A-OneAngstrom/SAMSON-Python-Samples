import datetime
import moviepy.editor as mpy

def create_gif(filenames, fps):
        '''
        function that creates a gif-file from images thanks to moviepy package
        '''
        output_gif_file = 'moviepy-%sam.gif' % datetime.datetime.now().strftime('%Y-%M-%d-%H-%M-%S')  # generate name for the gif-file
        output_mp4_file = 'moviepy-%sam.mp4' % datetime.datetime.now().strftime('%Y-%M-%d-%H-%M-%S')  # generate name for the mp4-file
        clip = mpy.ImageSequenceClip(filenames, fps=fps)                    # create a clip
        clip.write_gif(output_gif_file)                                     # save clip to the gif-file
        clip.write_videofile(output_mp4_file)                               # save clip to the mp4-file

# Load the structure file
mol new "poscar.pdb"

# Set display settings
display projection orthographic
display rendermode GLSL
axes location off

# Set background to white
color Display Background 8

# Apply CPK representation and color atoms by element
mol modstyle 0 0 CPK
mol modcolor 0 0 Element ;# Default: color atoms by their element type

# Double the radius of gold atoms
set gold_sel [atomselect top "name Au"] ;# Select all gold atoms (AU is the common symbol)
$gold_sel set radius 5.0                ;# Initialize radius to 1.0
$gold_sel delete

# Adjust the initial orientation for a side view (XZ plane)
rotate x by -90 ;# Rotate the camera to view the cell from the side

# Parameters for rotation
set num_frames 360         ;# Total number of frames for a full 360-degree rotation
set angle_increment [expr 360.0 / $num_frames]
set output_prefix "frame"  ;# Prefix for output frames

# Rotate the camera around the Z-axis
for {set i 0} {$i < $num_frames} {incr i} {
    rotate y by $angle_increment
    render snapshot "${output_prefix}_[format "%04d" $i].tga"
}

# Completion message
puts "Rendering completed. Use ffmpeg to create the movie:"
puts "ffmpeg -framerate 30 -i ${output_prefix}_%04d.tga -c:v libx264 -pix_fmt yuv420p rotation_movie.mp4"


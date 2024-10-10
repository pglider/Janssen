gcc main.c -o out.exe -O3 -Wall -lm
.\out.exe
magick mogrify -format jpg -background white -alpha remove -quality 95 OUTPUT/p*.eps
ffmpeg -framerate 15 -y -i C:\Users\Utilisateur\Documents\Janssen\OUTPUT\p%%05d.jpg -c:v libx264 C:\Users\Utilisateur\Documents\Janssen\OUTPUT\movie.avi
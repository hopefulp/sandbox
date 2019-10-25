### global

### initialization
set filename [molinfo top get filename]
if { $filename == "" } {
	puts "No loaded molecule."
	exit 1
}
#puts "** the filename of the loaded molecule: $filename"

if { [catch { exec grep CRYS $filename } msg] } {
	puts "No PBC information on $filename"
}

set a [join $msg " "]
set b [lrange $a 1 6]

pbc set "{ $b }"
puts "got PBC information from the BGF file: PBC set to { $b }"

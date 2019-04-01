use 5.14.0;

my @stuff = qw(just a list);

my %h = (one => \@stuff, two => "beta");

print(join(" ", keys(%h)), "\n");

my $three = $h{three};
my $sr = $h{one};

if($three) {
say("\$three is true");
}
if($sr) {
say("\$sr is true");
say(ref($sr));
}

print(join(" ", keys(%h)), "\n");

exit;


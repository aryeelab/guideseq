import svgwrite
import sys

boxWidth = 10
box_size = 15
v_spacing = 3

colors = {'G': '#F5F500', 'A': '#FF5454', 'T': '#00D118', 'C': '#26A8FF', 'N': '#B3B3B3'}

def parseSitesFile(infile):
    offtargets = []
    with open(infile, 'r') as f:
        f.readline()
        for line in f:
            line_items = line.split('\t')
            offtarget_sequence = line_items[21]
            offtarget_reads = line_items[15]
            ref_seq = line_items[32]
            if offtarget_sequence != '':
                offtargets.append({'seq': offtarget_sequence.strip(),
                                   'reads': int(offtarget_reads.strip())})
    offtargets = sorted(offtargets, key=lambda x: x['reads'], reverse=True)
    return offtargets, ref_seq


def visualizeOfftargets(infile, outfile):

    # Get offtargets array from file
    offtargets, ref_seq = parseSitesFile(infile)

    # Initiate canvas
    dwg = svgwrite.Drawing(outfile + '.svg', profile='full')

    # Define top and left margins
    x_offset = 20
    y_offset = 30

    # Draw ticks
    tick_locations = [1, len(ref_seq)]
    for x in tick_locations:
        dwg.add(dwg.text(str(x), insert=(x_offset + (x - 1) * box_size + 2, y_offset - 2), style="font-size:10px; font-family:Courier"))

    # Draw reference sequence row
    for i, c in enumerate(ref_seq):
        y = y_offset
        x = x_offset + i * box_size
        dwg.add(dwg.rect((x, y), (box_size, box_size), fill=colors[c]))
        dwg.add(dwg.text(c, insert=(x + 3, y + box_size - 3), fill='black', style ="font-size:15px; font-family:Courier"))


    dwg.add(dwg.text('Reads', insert=(x_offset + box_size * len(ref_seq) + 16, y_offset + box_size - 3), style="font-size:15px; font-family:Courier"))

    # Draw aligned sequence rows
    y_offset += 10 # leave some extra space after the reference row
    for j, seq in enumerate(offtargets):
        y = y_offset + j * box_size
        for i, c in enumerate(seq['seq']):
            x = x_offset + i * box_size
            if c == ref_seq[i] or ref_seq[i] == 'N':
                dwg.add(dwg.text(u"\u2022", insert=(x + 3, 2 * box_size + y - 3), fill='black', style ="font-size:15px; font-family:Courier"))
            else:
                dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), fill=colors[c]))
                dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style ="font-size:15px; font-family:Courier"))

        dwg.add(dwg.text(str(seq['reads']), insert=(box_size * (len(ref_seq) + 1) + 20, y_offset + box_size * (j + 2) - 2), fill='black', style="font-size:15px; font-family:Courier"))

    dwg.save()


def main():
    if len(sys.argv) == 3:
        visualizeOfftargets(sys.argv[1], sys.argv[2])
    else:
        print 'Usage: python visualization.py INFILE OUTFILE'


if __name__ == '__main__':
    main()
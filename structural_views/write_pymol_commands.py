#Make text file with pymol commands

def main():
    outFile = 'pymol_commands.txt'
    crowe_mabs = ['2050', '2082', '2094', '2096', '2130', '2165', '2196', '2479', '2499', '2677', '2832']
    antibodies = ['CR3022']

    # which set of data to show, and how to color surfaces.
    metric = 'max' # `max` looks way better than `total`
    color_min = 'white'
    color_max = 'red'

    # dictionary of views here:
    views = {'view1': """\nset_view (\
0.179795116,   -0.585253060,    0.790621936,\
-0.975049794,    0.000042952,    0.221752375,\
-0.129811451,   -0.810812414,   -0.570666671,\
0.002255535,    0.000589155, -237.384719849,\
-32.869243622,   26.272405624,   17.597675323,\
-55707.390625000, 56183.121093750,  -20.000000000 )\n""",

    }

    # save png files of views?
    save_png = False

    # make full list of antibodies
    for mab in crowe_mabs:
        ab = 'COV2-'+mab
        antibodies.append('COV2-'+mab)

    # begin writing to output file
    f = open(outFile, "w")
    f.write('#commands to load pdbs for antibody mapping\n\n')

    # set sequence view to off
    f.write('set seq_view, 0\n\n')

    # load pdb file for every antibody
    for ab in antibodies:
        f.write(f'load ../results/pdb_outputs/{ab}_400_6m0j_{metric}_escape.pdb\n')

    # rename structure for every antibody
    f.write('\n')
    for ab in antibodies:
        f.write(f'set_name {ab}_400_6m0j_{metric}_escape, {ab}_{metric}\n')

    # one-time commands:
    f.write('\nhide all\n')
    f.write(f'create ACE2, {antibodies[0]}_{metric} and chain A\n')

    # remove chain A (ACE2) from each structure.
    for ab in antibodies:
        f.write(f'remove {ab}_{metric} and chain A\n')

    # show ACE2 as gray cartoon
    f.write('\nshow cartoon, ACE2; color gray20, ACE2; set cartoon_transparency, 0.5, ACE2\n\n')

    # color surface representation by b-factor (which has been recoded to {metric}_escape)
    for ab in antibodies:
        f.write(f'show surface, {ab}_{metric}; show sticks, {ab}_{metric} and resn NAG; spectrum b, {color_min} {color_max}, {ab}_{metric}, minimum=0\n')

    # show each as surface
    f.write('\nhide all\n')

    for ab in antibodies:
        f.write(f'show surface, {ab}_{metric}\n')

    # change coloring of RBD_bind and RBD_express to blue white (where blue is negative)
    f.write(f'\nshow surface, RBD_bind; spectrum b, blue white, RBD_bind, minimum=-2, maximum=0; show sticks, RBD_bind and resn NAG')
    f.write(f'\nshow surface, RBD_expr; spectrum b, blue white, RBD_expr, minimum=-2, maximum=0; show sticks, RBD_expr and resn NAG\n')

    # iterate through each view we define above, get view, show as surface, take picture
    for view_n in views.keys():

        f.write(f'{views[view_n]}\n') # set one view you like

        # now get pictures of each
        for ab in antibodies:
            if save_png:
                f.write('\nhide all\n')
            f.write(f'show surface, {ab}_{metric}\n')

            if save_png:
                f.write(f'png ../structural_views/{ab}_{metric}_{view_n}.png, ray=1, 600, 600\n')
        if not save_png:
            f.write('\nsave surface_escape.pse')
            break

    f.close()

main()

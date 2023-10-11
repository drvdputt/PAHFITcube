#!/usr/bin/env python3
from pahfit.scripts.run_pahfit import initialize_parser
from pahfitcube.cube_model import CubeModel
from pahfit.model import Model
from pahfitcube.iohacks import read_cube
from pathlib import Path


def run_pahfit_cube(args):
    """Main function, called after parsing the args.

    Can also be run by importing it, and passing the right args. Run
    ./run_pahfit_cube.py --help to see which arguments can be specified.
    The ones currently required are 'spectrumfile', 'j'

    """
    # load spectral cube
    spec, _, wcs = read_cube(args.spectrumfile)
    spec.meta["instrument"] = args.instrumentname

    # set up cube model with single PAHFIT model as starting point
    model = Model.from_yaml(args.packfile)
    cube_model = CubeModel(model)

    # output directory
    output_dir = Path(args.o).resolve()
    output_dir.mkdir(exist_ok=True)
    # output file prefix (removes the extension .fits or .pickle)
    output_base = Path(args.spectrumfile).name.rsplit(".", 1)[0]

    # determine location of checkpoint files (if checkpointing is
    # desired)
    if args.resume:
        checkpoint_dir = output_dir / "checkpoints"
        checkpoint_dir.mkdir(exist_ok=True)
        checkpoint_prefix = str(checkpoint_dir / output_base)
    else:
        checkpoint_prefix = None

    # do the fit
    cube_model.fit(spec, checkpoint_prefix, maxiter=args.fit_maxiter)

    # save result
    cube_model.maps.save(wcs, str(output_dir / (output_base + "_maps.fits")))

    # # run everything
    # num_fits = len(spaxels)
    # if args.j > 1:
    #     with Pool(args.j) as p:
    #         parallel_iterator = p.imap(
    #             fit_spaxel_wrapper,
    #             ((spaxels[i], args) for i in range(num_fits)),
    #         )
    #         models = []
    #         for i, model in enumerate(parallel_iterator):
    #             models.append(model)
    #             print(f"Finished fit {i}/{num_fits}")
    # else:
    #     models = [fit_spaxel(s, args) for s in spaxels]


def main(args_list=None):
    """
    args_list: list of str
        override the command line arguments (useful for calling main
        from another script)
    """
    parser = initialize_parser()
    # Need to figure out the common parts between this set of arguments
    # and those set in initialize_parser, so that the single spectrum
    # and the cube case can each have there own arguments, without
    # unnecessary duplication.
    parser.add_argument("-j", type=int, default=1, help="Number of parallel processes")
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Load pmodel for pixel from file if present",
    )
    parser.add_argument("-o", default=".", help="Output directory")
    if args_list is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args_list)

    run_pahfit_cube(args)


if __name__ == "__main__":
    main()

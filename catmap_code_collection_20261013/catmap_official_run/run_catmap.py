"""Run CatMAP 0.3.1 input parsing and static endpoint-barrier checks only.

This entry point deliberately does not call ReactionModel.run(), so it does not
perform steady-state, coverage, rate, TOF, or mapper calculations.
"""

from catmap import ReactionModel, __version__

SURFACES = [
    'tife', 'timn', 'tihf', 'tire', 'timo', 'tiv', 'tizr',
    'tico', 'titi', 'tiw', 'tita', 'ticr', 'ti', 'tini',
]


def main():
    if __version__ != '0.3.1':
        raise RuntimeError('Expected local CatMAP 0.3.1; found ' + __version__)

    print('CatMAP', __version__)
    for surface in SURFACES:
        # Construction invokes TableParser and CatMAP's official verify().
        model = ReactionModel(setup_file='alkene_epoxidation.mkm')
        model.surface_names = [surface]
        model.parse()
        model.verify()
        energy_dict = model.scaler.get_electronic_energies([333.15, 1.0])
        energy_dict.setdefault('s', 0.0)
        if len(model.elementary_rxns) != 6:
            raise RuntimeError('Expected the six approved reaction expressions')
        print('surface={}; TableParser+atom/site balance: PASS'.format(surface))

        if surface == 'tire':
            # CatMAP's official static endpoint-barrier calculation; no rates.
            for index, rxn in enumerate(model.elementary_rxns, start=1):
                delta_g, forward_barrier = model.get_rxn_energy(rxn, energy_dict)
                _, reverse_barrier = model.get_rxn_energy(
                    model.reverse_rxn(rxn), energy_dict
                )
                print(
                    'TiRe step={}; dG_eV={:.9f}; Eforward_eV={:.9f}; '
                    'Ereverse_eV={:.9f}'.format(
                        index, delta_g, forward_barrier, reverse_barrier
                    )
                )

    print('ReactionModel.run(): NOT CALLED')


if __name__ == '__main__':
    main()

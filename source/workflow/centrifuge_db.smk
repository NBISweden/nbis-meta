centrifuge_db_input = expand(opj(config["centrifuge_dir"],
                                 "{base}.{i}.cf"),
                             i=[1,2,3], base=config["centrifuge_base"])

include: "../rules/Classify/prepare_centrifuge.smk"
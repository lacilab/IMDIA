# Generate command-line flags from argument dictionary
def gen_flags(args):
    flags = []

    for name, value in args.items():
        if value != False and value is not None:
            flag = "--" + name
            flags.append(flag)
            if value != True:
                flags.append(str(value))

    return flags
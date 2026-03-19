import sys
import argparse
from validation.validation_report import run_validation


def main():
    parser = argparse.ArgumentParser(description="Cryptic IP Binding Sites Pipeline")
    parser.add_argument(
        "--mode", default="validation", choices=["validation"], help="Mode to run in"
    )
    args = parser.parse_args()

    if args.mode == "validation":
        run_validation()
    else:
        print(f"Unknown mode: {args.mode}")
        sys.exit(1)


if __name__ == "__main__":
    main()

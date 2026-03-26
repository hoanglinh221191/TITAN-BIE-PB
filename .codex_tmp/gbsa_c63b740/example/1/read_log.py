import sys
import codecs

def main():
    try:
        with codecs.open('out.log', 'r', 'utf-16le') as f:
            lines = f.readlines()
            for line in lines:
                print(line.strip())
    except Exception as e:
        print(f"Error: {e}")

if __name__ == '__main__':
    main()

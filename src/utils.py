import pickle

def save_dataframe(df, filename):
    with open(filename, 'wb') as file:
        pickle.dump(df, file)


def load_dataframe(filename):
    with open(filename, 'rb') as file:
        return pickle.load(file)

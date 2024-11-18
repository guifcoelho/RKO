import rkopy

if __name__ == '__main__':

    def decoder(keys: list[float]):
        # print(keys)
        return sum(keys)

    best_sol = rkopy.solve(
        "Teste",
        1,
        False,
        "ParametersOffline.txt",
        True,
        10,
        decoder
    )

    print(best_sol)
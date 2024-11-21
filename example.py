import rkopy

if __name__ == '__main__':

    def decoder(keys: list[float]):
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

    print("\nBest solution:")
    print(best_sol)

    print(f"\nBest value: {decoder(best_sol)}")

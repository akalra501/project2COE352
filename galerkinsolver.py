import numpy as np
import matplotlib.pyplot as plt

def computeExactSolution(x, time):
    """
    Compute the exact solution for the problem.

    """
    return np.exp(-time) * np.sin(np.pi * x)

def computeSourceTerm(x, time):
    """
    Compute the source term for the problem.

    """
    return (np.pi**2 - 1) * np.exp(-time) * np.sin(np.pi * x)


def computeInitialCondition(x):
    """
    Compute the initial condition for the problem.

    """
    return np.sin(np.pi * x)

def createElementMapping(numNodes):
    """
    Create element-to-node mapping for FEM.

    """
    return np.vstack((np.arange(0, numNodes - 1), np.arange(1, numNodes))).T

def createEmptyMatrices(numNodes, numTimesteps):
    """
    Create empty matrices for stiffness, mass, and force.

    """
    return (np.zeros((numNodes, numNodes)),
            np.zeros((numNodes, numNodes)),
            np.zeros((numNodes, numTimesteps + 1)))

def createBasisFunctions(elementLength):
    """
    Create basis functions and associated data for FEM.

    """
    basisFunc1 = lambda zeta: (1 - zeta) / 2
    basisFunc2 = lambda zeta: (1 + zeta) / 2
    
    basisFunctionDerivatives = np.array([-1 / 2, 1 / 2])
    quadraturePoints = [-1 / np.sqrt(3), 1 / np.sqrt(3)]
    basisFunctionsAtQuad = np.array([[basisFunc1(zeta), basisFunc2(zeta)] for zeta in quadraturePoints])
    return (basisFunctionsAtQuad, basisFunctionDerivatives, 2 / elementLength, elementLength / 2)

def assembleFEMMatrices(numNodes, numTimesteps, stiffnessMatrix, massMatrix, forceMatrix, mapping, basisFunctions, 
                        derivatives, derivativeScaling, integralScaling, elementLength):
    """
    Assemble FEM matrices for the problem.

    """
    quadraturePoints = [-0.57735026919, 0.57735026919]
    for elementIndex in range(numNodes - 1):
        localMassMatrix = np.zeros((2, 2))
        localStiffnessMatrix = np.zeros((2, 2))
        for i in range(2):
            for j in range(2):
                localMassMatrix[i, j] = sum(basisFunctions[i, k] * basisFunctions[j, k] for k in range(2)) * elementLength 
                localStiffnessMatrix[i, j] = derivatives[i] * derivativeScaling * derivatives[j] * derivativeScaling*integralScaling * 2
        globalNodes = mapping[elementIndex].astype(int)
        for i in range(2):
            for j in range(2):
                stiffnessMatrix[globalNodes[i], globalNodes[j]] += localStiffnessMatrix[i, j]
                massMatrix[globalNodes[i], globalNodes[j]] += localMassMatrix[i, j]
        forceMatrix[elementIndex, :] = -sum(computeSourceTerm(zeta, numTimesteps) * basisFunctions[0, k] for k, 
                                            zeta in enumerate(quadraturePoints)) * (1 / 8)
    return massMatrix, stiffnessMatrix, forceMatrix

def applyDirichletConditions(massMatrix, numNodes):
    """
    Apply Dirichlet boundary conditions to the mass matrix.

    """
    massMatrix[0, :] = massMatrix[-1, :] = massMatrix[:, 0] = massMatrix[:, -1] = 0
    massMatrix[0, 0] = massMatrix[-1, -1] = 1
    dirichletBoundaryConditions = np.eye(numNodes)
    dirichletBoundaryConditions[0, 0] = dirichletBoundaryConditions[-1, -1] = 0
    return massMatrix, dirichletBoundaryConditions

# Time integration functions
def setupEulerMatrices(massMatrix, stiffnessMatrix, timestepSize):
    """
    Set up matrices for Euler time integration methods.

    """
    inverseMassMatrix = np.linalg.inv(massMatrix)
    eulerMatrix = (1 / timestepSize) * massMatrix + stiffnessMatrix
    return inverseMassMatrix, np.dot(inverseMassMatrix, stiffnessMatrix), np.linalg.inv(eulerMatrix)

def solveEulerTimesteps(numNodes, numTimesteps, timestepSize, massStiffnessProduct, inverseMassMatrix, massMatrix, forceMatrix, boundaryConditions, eulerMethod, nodeCoords, inverseEulerMatrix):
    """
    Solve the problem using Euler time integration methods.

    """
    solution = np.zeros((numNodes, numTimesteps + 1))
    solution[:, 0] = computeInitialCondition(nodeCoords)
    for t in range(numTimesteps):
        if eulerMethod == 'FE':
            solution[:, t + 1] = solution[:, t] - timestepSize * massStiffnessProduct.dot(solution[:, t]) + timestepSize * inverseMassMatrix.dot(forceMatrix[:, t])
        elif eulerMethod == 'BE':
            solution[:, t + 1] = (1 / timestepSize) * inverseEulerMatrix.dot(massMatrix.dot(solution[:, t])) + inverseEulerMatrix.dot(forceMatrix[:, t])
        solution[:, t + 1] = boundaryConditions.dot(solution[:, t + 1])
    return solution

# Visualization
def plotComparison(continuousCoords, analyticalSolution, femCoords, numericalSolution, numTimesteps, eulerMethod):
    """
    Plot comparison between analytical and numerical solutions.

    """
    plt.plot(continuousCoords, analyticalSolution, label='Analytical Solution', color="black")
    methodLabel = "Forward Euler" if eulerMethod == "FE" else "Backward Euler"
    plt.plot(femCoords, numericalSolution[:, numTimesteps], label=f'{methodLabel} (t = {numTimesteps})', color="red")
    plt.xlabel('x')
    plt.ylabel('Solution')
    plt.title('Plot of Analytical vs Numerical Solution')
    plt.legend()
    plt.show()

def solveFEMWithMethod(numNodes, numTimesteps, nodeCoords, timeSteps, elementLength, timestepSize, method):
    """
    Solve the FEM problem for the specified method.

    """
    massMatrix, stiffnessMatrix, forceMatrix = createEmptyMatrices(numNodes, numTimesteps)
    mapping = createElementMapping(numNodes)
    basisFunctions, derivatives, derivativeScaling, integralScaling = createBasisFunctions(elementLength)
    massMatrix, stiffnessMatrix, forceMatrix = assembleFEMMatrices(
        numNodes, timeSteps, stiffnessMatrix, massMatrix, forceMatrix,
        mapping, basisFunctions, derivatives, derivativeScaling, integralScaling, elementLength
    )
    massMatrix, dirichletBC = applyDirichletConditions(massMatrix, numNodes)
    inverseMassMatrix, massStiffnessProduct, inverseEulerMatrix = setupEulerMatrices(massMatrix, stiffnessMatrix, timestepSize)
    solution = solveEulerTimesteps(
        numNodes, numTimesteps, timestepSize, massStiffnessProduct,
        inverseMassMatrix, massMatrix, forceMatrix, dirichletBC, method,
        nodeCoords, inverseEulerMatrix
    )
    continuousCoords = np.linspace(0, 1, 500)
    analyticalSolution = computeExactSolution(continuousCoords, 1)
    plotComparison(continuousCoords, analyticalSolution, nodeCoords, solution, numTimesteps, method)


def chooseEulerMethod():
    """
    Prompt the user to choose an Euler method and ensure valid input.

    """
    while True:
        method = input("Type FE for Forward Euler and BE for Backward Euler (NO OTHER VALUES): ").upper()
        if method in ['FE', 'BE']:
            return method
        print("Invalid choice. Please type FE or BE.")

def main():
    """
    Main function to solve the FEM problem.

    """
    try:
        numNodes = int(input("Input the number of spatial nodes(N): "))
        numTimesteps = int(input("Input the number of timesteps(t): "))
        nodeCoords = np.linspace(0, 1, numNodes)
        elementLength = nodeCoords[1] - nodeCoords[0]
        timestepSize = 1 / numTimesteps
        timeSteps = np.linspace(0, 1, numTimesteps + 1)
        method = chooseEulerMethod()
        solveFEMWithMethod(numNodes, numTimesteps, nodeCoords, timeSteps, elementLength, timestepSize, method)
    except ValueError:
        print("This input is invalid, please enter only numeric values.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()

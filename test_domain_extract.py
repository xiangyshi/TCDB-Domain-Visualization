import subprocess
import sys


def run_test(args, description):
    print(f"\n{'='*30}\nRunning: {description}\n{'='*30}")
    result = subprocess.run(
        [sys.executable, 'domain_extract.py'] + args,
        capture_output=True,
        text=True
    )
    print("STDOUT:\n", result.stdout)
    print("STDERR:\n", result.stderr)
    print(f"Exit code: {result.returncode}\n")
    return {
        'description': description,
        'args': args,
        'exit_code': result.returncode,
        'stdout': result.stdout,
        'stderr': result.stderr
    }


def main():
    test_results = []
    
    # Test 1: CDD input
    test_results.append(run_test(['-c', 'singles.cdd', '-t', '1.A.12'], 
                                "CDD input: -c singles.cdd -t 1.A.12"))

    # Test 2: Rescue input
    test_results.append(run_test(['-r', 'rescued', '-t', '1.A.12'], 
                                "Rescue input: -r rescued -t 1.A.12"))

    # Print summary
    print("\n" + "="*50)
    print("TEST SUMMARY")
    print("="*50)
    
    successful_tests = [test for test in test_results if test['exit_code'] == 0]
    failed_tests = [test for test in test_results if test['exit_code'] != 0]
    
    print(f"\nSuccessful Tests ({len(successful_tests)}):")
    for test in successful_tests:
        print(f"✓ {test['description']}")
    
    print(f"\nFailed Tests ({len(failed_tests)}):")
    for test in failed_tests:
        print(f"✗ {test['description']} (Exit code: {test['exit_code']})")
    
    print(f"\nTotal Tests: {len(test_results)}")
    print(f"Passed: {len(successful_tests)}")
    print(f"Failed: {len(failed_tests)}")


if __name__ == "__main__":
    main()
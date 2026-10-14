"""Run and audit the appendix outside its published core code."""
from pathlib import Path
import ast, contextlib, hashlib, importlib.util, io, json, math
ROOT = Path(__file__).resolve().parent


def check():
    path = ROOT / 'appendices/solve.py'
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    target = ROOT / 'build/revision/appendix_validation.json'
    if target.exists():
        cached = json.loads(target.read_text(encoding='utf-8'))
        if cached.get('source_sha256') == digest:
            return cached
    spec = importlib.util.spec_from_file_location('appendix_solver', path)
    core = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(core)
    area = core.calibrate()
    # Execute the actual entry point, reusing the just-computed calibration.
    core.calibrate = lambda: area
    capture = io.StringIO()
    with contextlib.redirect_stdout(capture):
        computed = core.main()
    lines = capture.getvalue().strip().splitlines()
    assert len(lines) == 3
    heights = [float(line) * 1000 for line in lines]
    r = json.loads((ROOT / 'results.json').read_text(encoding='utf-8'))
    expected = [r['free_g_extrapolation']['maximum_height_m'],
                r['maximum_height_m'], r['spin_response']['maximum_height_m']]
    differences = [abs(a-b) for a,b in zip(heights, expected)]
    assert max(differences) < 3., differences
    assert [round(h/1000,2) for h in heights] == [round(h/1000,2) for h in expected]
    assert abs(area-r['calibration']['area_m2']) < .001
    residuals = [core.maximum(h, area, k)-limit for h,k,limit in
                 zip(heights, ['load','heat','spin'], [5.,390.,5.])]
    assert max(map(abs,residuals)) < .001
    try:
        core.atmosphere(87000)
        raise AssertionError('Default atmosphere did not reject extrapolation')
    except ValueError:
        pass
    for forbidden in ['matplotlib','sensitivity','rmse','unittest']:
        assert forbidden not in path.read_text(encoding='utf-8-sig').lower()
    tree = ast.parse(path.read_text(encoding='utf-8'))
    imports = [n.module if isinstance(n,ast.ImportFrom) else n.names[0].name
               for n in ast.walk(tree) if isinstance(n,(ast.Import,ast.ImportFrom))]
    assert sorted(imports) == ['bisect','math']
    vacuum = core.integrate(lambda t,y:(-y[1],core.G0),(1000.,0.),0.)
    assert abs(vacuum.times[-1]-math.sqrt(2000/core.G0)) < 1e-7
    assert abs(vacuum.states[-1][1]-math.sqrt(2000*core.G0)) < 1e-6
    constant = core.integrate(lambda t,y:(-y[1],core.G0-.01*y[1]**2),(1000.,0.),0.)
    terminal = math.sqrt(core.G0/.01)
    analytic_speed = terminal*math.tanh(core.G0*constant.times[-1]/terminal)
    assert abs(constant.states[-1][1]-analytic_speed)<1e-6
    descent = computed['descent']; chute = descent['parachute']
    assert chute['trajectory'].states[0] == descent['freefall'].states[-1]
    assert chute['trajectory'].states[-1][0] == 0 and chute['inflation_completed']
    assert 0 < chute['full_inflation_height_m'] < core.OPENING
    import model
    reference = model.chute(core.OPENING,chute['opening_speed_mps'],core.MASS,area)
    comparison = {'opening_load_g':abs(chute['opening_load_g']-reference['summary']['opening_g_peak']),
                  'landing_speed_mps':abs(chute['landing_speed_mps']-reference['summary']['landing_mps'])}
    assert max(comparison.values()) < 1e-5,comparison
    assert chute['opening_load_g'] < 5 and chute['landing_speed_mps'] < 6
    result = {'status':'PASS','source_sha256':digest,'area_m2':area,
              'stdout':capture.getvalue(),'heights_m':heights,
              'height_difference_m':differences,'boundary_residuals':residuals,
              'standalone':True,'three_outputs_only':True,'imports':imports,
              'internal_adaptive_rk4_verified':True,
              'canopy_reference_differences':comparison,
              'canopy':{k:v for k,v in chute.items() if k!='trajectory'},
              'total_time_s':descent['total_time_s'],
              'opening_handoff_and_ground_event_verified':True}
    target.write_text(json.dumps(result,indent=2),encoding='utf-8')
    (ROOT/'results/appendix_three_heights.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    return result


if __name__ == '__main__':
    print(json.dumps(check(),indent=2))
